use proptest::prelude::*;

use crate::elliptic_curves::{
    AffineCurveModel, AffinePoint, CurveError, GroupCurveModel, ShortWeierstrassCurve,
    ShortWeierstrassFunction, ShortWeierstrassFunctionField, ShortWeierstrassFunctionFieldPoint,
};
use crate::fields::{AmbientField, Field, Fp, PthRootExtraction, Q, RationalFunction};
use crate::polynomials::DensePolynomial;
use crate::polynomials::evaluation::evaluate_dense;
use crate::proptest_support::config::{CurveStrategyConfig, PolynomialStrategyConfig};
use crate::proptest_support::elliptic_curves::{
    arb_short_weierstrass_function_case, arb_short_weierstrass_function_pair_case,
};

type F17 = Fp<17>;

fn f17_dense(values: &[u64]) -> DensePolynomial<F17> {
    DensePolynomial::<F17>::new(values.iter().copied().map(F17::elem_from_u64).collect())
}

fn q_dense(values: &[(i64, i64)]) -> DensePolynomial<Q> {
    DensePolynomial::<Q>::new(
        values
            .iter()
            .map(|&(numerator, denominator)| {
                Q::div(&Q::from_i64(numerator), &Q::from_i64(denominator))
                    .expect("denominator should be non-zero")
            })
            .collect(),
    )
}

fn f17_curve() -> ShortWeierstrassCurve<F17> {
    ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(2), F17::elem_from_u64(3))
        .expect("curve should be nonsingular")
}

fn q_curve() -> ShortWeierstrassCurve<Q> {
    ShortWeierstrassCurve::<Q>::new(Q::from_i64(-1), Q::from_i64(0)).expect("curve should exist")
}

#[test]
fn function_field_family_exposes_zero_one_x_and_y() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());

    assert!(field.zero().is_zero());
    assert!(field.one().is_one());
    assert_eq!(
        field.x().a_part(),
        &RationalFunction::<F17>::indeterminate()
    );
    assert!(field.x().b_part().is_zero());
    assert!(field.y().a_part().is_zero());
    assert!(field.y().b_part().is_one());
}

#[test]
fn generic_point_is_the_pair_of_distinguished_coordinate_functions() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(generic.x(), Some(&field.x()));
    assert_eq!(generic.y(), Some(&field.y()));
}

#[test]
fn translating_the_generic_point_by_infinity_is_the_identity() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let translated = field
        .translate_generic_point_by_finite_point(&AffinePoint::Infinity)
        .expect("identity translation should succeed");

    assert_eq!(translated, field.generic_point());
}

#[test]
fn adding_infinity_is_neutral_for_function_field_points() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(
        field.add_points(&generic, &ShortWeierstrassFunctionFieldPoint::Infinity),
        Ok(generic.clone())
    );
    assert_eq!(
        field.add_points(&ShortWeierstrassFunctionFieldPoint::Infinity, &generic),
        Ok(generic)
    );
}

#[test]
fn adding_a_point_and_its_negative_gives_infinity() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let point = field.generic_point();
    let negated = field.neg_point(&point).expect("negation should succeed");

    assert_eq!(
        field.add_points(&point, &negated),
        Ok(ShortWeierstrassFunctionFieldPoint::Infinity)
    );
}

#[test]
fn adding_a_constant_point_to_the_generic_point_matches_translation_helper() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let translation_point = curve
        .point(Fp::<41>::from_i64(40), Fp::<41>::from_i64(0))
        .expect("sample translation point should lie on the curve");
    let generic = field.generic_point();
    let constant = field
        .from_affine_point(&translation_point)
        .expect("constant embedding should succeed");

    assert_eq!(
        field.add_points(&generic, &constant),
        field.translate_generic_point_by_finite_point(&translation_point)
    );
}

#[test]
fn doubling_the_generic_point_matches_the_generic_doubling_helper() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let generic = field.generic_point();

    assert_eq!(field.double_point(&generic), field.double_generic_point());
}

#[test]
fn generic_scalar_multiplication_by_three_matches_curve_scalar_multiplication_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let expected = curve
        .mul_scalar(&sample_point, 3)
        .expect("scalar multiplication should stay on the curve");
    let triple = field
        .generic_point_multiple(3)
        .expect("generic scalar multiplication should succeed");
    let triple_x = triple.x().expect("triple stays affine");
    let triple_y = triple.y().expect("triple stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(triple_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(triple_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}

#[test]
fn translating_the_generic_point_matches_curve_addition_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let translation_point = curve
        .point(Fp::<41>::from_i64(40), Fp::<41>::from_i64(0))
        .expect("sample translation point should lie on the curve");
    let expected = curve
        .add(&sample_point, &translation_point)
        .expect("sample translation should stay on the curve");
    let translated = field
        .translate_generic_point_by_finite_point(&translation_point)
        .expect("generic translation should succeed");
    let translated_x = translated.x().expect("translation stays affine");
    let translated_y = translated.y().expect("translation stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(&translated_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(&translated_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}

#[test]
fn doubling_the_generic_point_matches_curve_doubling_on_a_sample_point() {
    let curve = ShortWeierstrassCurve::<Fp<41>>::new(Fp::<41>::from_i64(2), Fp::<41>::from_i64(3))
        .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<Fp<41>>::new(curve.clone());
    let sample_point = curve
        .point(Fp::<41>::from_i64(3), Fp::<41>::from_i64(6))
        .expect("sample point should lie on the curve");
    let expected = curve
        .double(&sample_point)
        .expect("sample doubling should stay on the curve");
    let doubled = field
        .double_generic_point()
        .expect("generic doubling should succeed");
    let doubled_x = doubled.x().expect("doubling stays affine");
    let doubled_y = doubled.y().expect("doubling stays affine");

    assert_eq!(
        evaluate_short_weierstrass_function_at_point(&doubled_x, &sample_point),
        AffinePoint::x_coordinate(&expected).cloned()
    );
    assert_eq!(
        evaluate_short_weierstrass_function_at_point(&doubled_y, &sample_point),
        AffinePoint::y_coordinate(&expected).cloned()
    );
}

#[test]
fn function_value_embeds_rational_functions_in_the_a_part() {
    let curve = f17_curve();
    let rational = RationalFunction::<F17>::new(f17_dense(&[1]), f17_dense(&[1, 1]))
        .expect("rational function should exist");
    let function = ShortWeierstrassFunction::<F17>::from_rational_function(curve, rational.clone());

    assert_eq!(function.a_part(), &rational);
    assert!(function.b_part().is_zero());
}

#[test]
fn conjugate_flips_only_the_y_component() {
    let curve = f17_curve();
    let function = ShortWeierstrassFunction::<F17>::new(
        curve,
        RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 2])),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[3])),
    );

    let conjugate = function.conjugate();
    assert_eq!(conjugate.a_part(), function.a_part());
    assert_eq!(conjugate.b_part(), &function.b_part().neg());
}

#[test]
fn norm_matches_a_squared_minus_f_b_squared() {
    let curve = f17_curve();
    let function = ShortWeierstrassFunction::<F17>::new(
        curve.clone(),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[1, 1])),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
    );

    let expected = RationalFunction::<F17>::from_polynomial(f17_dense(&[15, 0, 1, 16]));
    assert_eq!(function.norm(), expected);
}

#[test]
fn multiplication_uses_the_short_weierstrass_relation() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let x = field.x();
    let y = field.y();

    let y_squared = y.mul(&y).expect("same-curve multiplication should work");
    let rhs = ShortWeierstrassFunction::<F17>::from_rational_function(
        field.curve().clone(),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1])),
    );

    assert_eq!(y_squared, rhs);
    assert_eq!(
        x.mul(&y).expect("same-curve multiplication should work"),
        y.mul(&x).expect("commutative")
    );
}

#[test]
fn inverse_uses_conjugate_over_norm() {
    let field = ShortWeierstrassFunctionField::<Q>::new(q_curve());
    let one = field.one();
    let y = field.y();
    let element = one.add(&y).expect("same-curve addition should work");

    let inverse = element.inverse().expect("1 + y should be invertible");
    let product = element.mul(&inverse).expect("multiplication should work");

    assert!(product.is_one());
}

#[test]
fn derivative_of_x_is_one() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let derivative = field.x().derivative();

    assert!(derivative.is_one());
}

#[test]
fn derivative_of_y_matches_f_prime_over_two_y_in_basis_form() {
    let field = ShortWeierstrassFunctionField::<Q>::new(q_curve());
    let derivative = field.y().derivative();

    assert!(derivative.a_part().is_zero());
    assert_eq!(
        derivative.b_part(),
        &RationalFunction::<Q>::new(
            q_dense(&[(-1, 1), (0, 1), (3, 1)]),
            q_dense(&[(0, 1), (-2, 1), (0, 1), (2, 1)])
        )
        .expect("rational function should exist")
    );
}

#[test]
fn derivative_satisfies_product_rule_on_small_example() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let x = field.x();
    let y = field.y();
    let product = x.mul(&y).expect("same-curve multiplication should work");
    let left = product.derivative();
    let right = x
        .derivative()
        .mul(&y)
        .expect("same-curve multiplication should work")
        .add(
            &x.mul(&y.derivative())
                .expect("same-curve multiplication should work"),
        )
        .expect("same-curve addition should work");

    assert_eq!(left, right);
}

#[test]
fn short_weierstrass_function_pth_root_recovers_the_y_generator_from_y_to_the_p() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let y = field.y();
    let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));
    let y_to_the_p = ShortWeierstrassFunction::<F17>::new(
        field.curve().clone(),
        RationalFunction::<F17>::constant(F17::zero()),
        {
            let mut result = RationalFunction::<F17>::constant(F17::one());
            for _ in 0..8 {
                result = result.mul(&rhs);
            }
            result
        },
    );

    assert_eq!(y_to_the_p.pth_root(), Some(y));
    assert!(y_to_the_p.has_pth_root());
}

fn evaluate_short_weierstrass_function_at_point<F: Field>(
    function: &ShortWeierstrassFunction<F>,
    point: &AffinePoint<F>,
) -> Option<F::Elem> {
    let x = AffinePoint::x_coordinate(point)?;
    let y = AffinePoint::y_coordinate(point)?;
    let a_value = evaluate_rational_function_at_x(function.a_part(), x)?;
    let b_value = evaluate_rational_function_at_x(function.b_part(), x)?;

    Some(F::add(&a_value, &F::mul(y, &b_value)))
}

fn evaluate_rational_function_at_x<F: Field>(
    function: &RationalFunction<F>,
    x: &F::Elem,
) -> Option<F::Elem> {
    let numerator = evaluate_dense(function.numerator(), x).ok()?;
    let denominator = evaluate_dense(function.denominator(), x).ok()?;

    if F::is_zero(&denominator) {
        None
    } else {
        F::div(&numerator, &denominator).ok()
    }
}

#[test]
fn short_weierstrass_function_pth_root_rejects_the_x_generator() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());

    assert_eq!(field.x().pth_root(), None);
    assert!(!field.x().has_pth_root());
}

#[test]
fn short_weierstrass_function_pth_root_handles_rational_a_part_examples() {
    let curve = f17_curve();
    let function = ShortWeierstrassFunction::<F17>::from_rational_function(
        curve.clone(),
        RationalFunction::<F17>::new(
            DensePolynomial::<F17>::new({
                let mut coefficients = vec![F17::zero(); 18];
                coefficients[17] = F17::one();
                coefficients
            }),
            DensePolynomial::<F17>::new({
                let mut coefficients = vec![F17::zero(); 18];
                coefficients[0] = F17::one();
                coefficients[17] = F17::one();
                coefficients
            }),
        )
        .expect("example rational function should exist"),
    );

    let expected = ShortWeierstrassFunction::<F17>::from_rational_function(
        curve,
        RationalFunction::<F17>::new(f17_dense(&[0, 1]), f17_dense(&[1, 1]))
            .expect("x / (1 + x) should exist"),
    );

    assert_eq!(function.pth_root(), Some(expected));
    assert!(function.has_pth_root());
}

#[test]
fn short_weierstrass_function_pth_root_recovers_a_mixed_example() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let root = field
        .x()
        .add(&field.y())
        .expect("same-curve addition should work");
    let rhs = RationalFunction::<F17>::from_polynomial(f17_dense(&[3, 2, 0, 1]));
    let function = ShortWeierstrassFunction::<F17>::new(
        field.curve().clone(),
        RationalFunction::<F17>::from_polynomial({
            let mut coefficients = vec![F17::zero(); 18];
            coefficients[17] = F17::one();
            DensePolynomial::<F17>::new(coefficients)
        }),
        {
            let mut result = RationalFunction::<F17>::constant(F17::one());
            for _ in 0..8 {
                result = result.mul(&rhs);
            }
            result
        },
    );

    assert_eq!(function.pth_root(), Some(root));
}

#[test]
fn evaluate_polynomial_in_x_substitutes_inside_the_function_field() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let x_plus_y = field
        .x()
        .add(&field.y())
        .expect("same-curve addition should work");
    let polynomial = f17_dense(&[3, 2, 1]);
    let expected = x_plus_y
        .mul(&x_plus_y)
        .expect("same-curve multiplication should work")
        .add(
            &x_plus_y
                .mul(&ShortWeierstrassFunction::<F17>::from_rational_function(
                    field.curve().clone(),
                    RationalFunction::<F17>::constant(F17::elem_from_u64(2)),
                ))
                .expect("same-curve multiplication should work"),
        )
        .expect("same-curve addition should work")
        .add(&ShortWeierstrassFunction::<F17>::from_rational_function(
            field.curve().clone(),
            RationalFunction::<F17>::constant(F17::elem_from_u64(3)),
        ))
        .expect("same-curve addition should work");

    assert_eq!(
        ShortWeierstrassFunction::<F17>::evaluate_polynomial_in_x(&polynomial, &x_plus_y)
            .expect("polynomial substitution should work"),
        expected
    );
}

#[test]
fn substitute_rational_function_in_x_embeds_regular_rational_substitution() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let function = RationalFunction::<F17>::new(f17_dense(&[1, 0, 1]), f17_dense(&[1, 1]))
        .expect("denominator should be non-zero");

    assert_eq!(
        ShortWeierstrassFunction::<F17>::substitute_rational_function_in_x(&function, &field.x())
            .expect("substitution at x should work"),
        field.from_rational_function(function)
    );
}

#[test]
fn substitute_rational_function_in_x_rejects_zero_denominator_after_substitution() {
    let curve = q_curve();
    let zero = ShortWeierstrassFunction::<Q>::zero(curve.clone());
    let function = RationalFunction::<Q>::new(q_dense(&[(1, 1)]), q_dense(&[(0, 1), (1, 1)]))
        .expect("denominator should be non-zero");

    assert_eq!(
        ShortWeierstrassFunction::<Q>::substitute_rational_function_in_x(&function, &zero),
        Err(CurveError::NonInvertibleFunctionFieldElement)
    );
}

#[test]
fn inverse_rejects_zero_norm_elements() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let zero = field.zero();

    assert_eq!(
        zero.inverse(),
        Err(CurveError::NonInvertibleFunctionFieldElement)
    );
}

#[test]
fn operations_reject_incompatible_curves() {
    let first_curve = f17_curve();
    let second_curve =
        ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(5), F17::elem_from_u64(7))
            .expect("curve should be nonsingular");

    let left = ShortWeierstrassFunction::<F17>::one(first_curve);
    let right = ShortWeierstrassFunction::<F17>::one(second_curve);

    assert_eq!(
        left.add(&right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
    assert_eq!(
        left.mul(&right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
    assert_eq!(
        left.div(&right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
}

#[test]
fn debug_output_mentions_curve_and_components() {
    let function = ShortWeierstrassFunction::<F17>::new(
        f17_curve(),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[1])),
        RationalFunction::<F17>::from_polynomial(f17_dense(&[0, 1])),
    );

    let debug = format!("{function:?}");
    assert!(debug.contains("ShortWeierstrassFunction"));
    assert!(debug.contains("curve_a"));
    assert!(debug.contains("a_part"));
    assert!(debug.contains("b_part"));
}

#[test]
fn ambient_field_zero_one_and_equality_match_the_function_field_family() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let zero = AmbientField::zero(&field);
    let one = AmbientField::one(&field);

    assert!(AmbientField::is_zero(&field, &zero));
    assert!(AmbientField::is_one(&field, &one));
    assert!(AmbientField::eq(&field, &field.x(), &field.x()));
}

#[test]
fn ambient_field_default_sub_and_div_work_for_function_fields() {
    let field = ShortWeierstrassFunctionField::<F17>::new(f17_curve());
    let x = field.x();
    let one = field.one();
    let x_plus_one = AmbientField::add(&field, &x, &one).expect("addition should work");
    let recovered_one =
        AmbientField::sub(&field, &x_plus_one, &x).expect("default subtraction should work");
    let recovered_x = AmbientField::div(
        &field,
        &x.mul(&one).expect("multiplication should work"),
        &one,
    )
    .expect("default division should work");

    assert!(AmbientField::eq(&field, &recovered_one, &one));
    assert!(AmbientField::eq(&field, &recovered_x, &x));
}

#[test]
fn ambient_field_reports_incompatible_curve_operations() {
    let first_curve = f17_curve();
    let second_curve =
        ShortWeierstrassCurve::<F17>::new(F17::elem_from_u64(5), F17::elem_from_u64(7))
            .expect("curve should be nonsingular");
    let field = ShortWeierstrassFunctionField::<F17>::new(first_curve.clone());
    let left = ShortWeierstrassFunction::<F17>::one(first_curve);
    let right = ShortWeierstrassFunction::<F17>::one(second_curve);

    assert_eq!(
        AmbientField::add(&field, &left, &right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
    assert_eq!(
        AmbientField::mul(&field, &left, &right),
        Err(CurveError::IncompatibleFunctionFieldCurves)
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn derivative_is_linear_over_same_curve_samples(
        case in arb_short_weierstrass_function_pair_case::<17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        let sum = case.left.add(&case.right).expect("generator should keep both functions on the same curve");
        let left = sum.derivative();
        let right = case.left
            .derivative()
            .add(&case.right.derivative())
            .expect("derivatives should stay on the same curve");

        prop_assert_eq!(left, right);
    }

    #[test]
    fn derivative_of_embedded_rational_constant_is_zero(
        case in arb_short_weierstrass_function_case::<17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        let constant = ShortWeierstrassFunction::<F17>::from_rational_function(
            case.curve.clone(),
            RationalFunction::<F17>::constant(F17::elem_from_u64(9)),
        );

        prop_assert!(constant.derivative().is_zero());
    }

    #[test]
    fn derivative_commutes_with_conjugation(
        case in arb_short_weierstrass_function_case::<17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        let left = case.function.derivative().conjugate();
        let right = case.function.conjugate().derivative();

        prop_assert_eq!(left, right);
    }
}
