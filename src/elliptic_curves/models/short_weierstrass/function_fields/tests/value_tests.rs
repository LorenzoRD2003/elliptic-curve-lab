use proptest::prelude::*;

use crate::elliptic_curves::{
    CurveError, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{Q, rational_function_field::RationalFunction};
use crate::proptest_support::{
    config::{CurveStrategyConfig, PolynomialStrategyConfig},
    elliptic_curves::{
        arb_short_weierstrass_function_case, arb_short_weierstrass_function_pair_case,
    },
};

use super::shared::{F17, f17_curve, f17_dense, q_curve, q_dense};

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
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
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
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            Q,
        >::new(q_curve());
    let one = field.one();
    let y = field.y();
    let element = one.add(&y).expect("same-curve addition should work");

    let inverse = element.inverse().expect("1 + y should be invertible");
    let product = element.mul(&inverse).expect("multiplication should work");

    assert!(product.is_one());
}

#[test]
fn evaluate_polynomial_at_function_x_substitutes_inside_the_function_field() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
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
                    RationalFunction::<F17>::constant(F17::from_i64(2)),
                ))
                .expect("same-curve multiplication should work"),
        )
        .expect("same-curve addition should work")
        .add(&ShortWeierstrassFunction::<F17>::from_rational_function(
            field.curve().clone(),
            RationalFunction::<F17>::constant(F17::from_i64(3)),
        ))
        .expect("same-curve addition should work");

    assert_eq!(
        ShortWeierstrassFunction::<F17>::evaluate_polynomial_at_function_x(&polynomial, &x_plus_y)
            .expect("polynomial substitution should work"),
        expected
    );
}

#[test]
fn substitute_rational_function_at_function_x_embeds_regular_rational_substitution() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let function = RationalFunction::<F17>::new(f17_dense(&[1, 0, 1]), f17_dense(&[1, 1]))
        .expect("denominator should be non-zero");

    assert_eq!(
        ShortWeierstrassFunction::<F17>::substitute_rational_function_at_function_x(
            &function,
            &field.x(),
        )
        .expect("substitution at x should work"),
        ShortWeierstrassFunction::<F17>::from_rational_function(f17_curve(), function)
    );
}

#[test]
fn substitute_rational_function_at_function_x_rejects_zero_denominator_after_substitution() {
    let curve = q_curve();
    let zero = ShortWeierstrassFunction::<Q>::zero(curve.clone());
    let function = RationalFunction::<Q>::new(q_dense(&[(1, 1)]), q_dense(&[(0, 1), (1, 1)]))
        .expect("denominator should be non-zero");

    assert_eq!(
        ShortWeierstrassFunction::<Q>::substitute_rational_function_at_function_x(&function, &zero),
        Err(CurveError::NonInvertibleFunctionFieldElement)
    );
}

#[test]
fn inverse_rejects_zero_norm_elements() {
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F17,
        >::new(f17_curve());
    let zero = field.zero();

    assert_eq!(
        zero.inverse(),
        Err(CurveError::NonInvertibleFunctionFieldElement)
    );
}

#[test]
fn operations_reject_incompatible_curves() {
    let first_curve = f17_curve();
    let second_curve = crate::elliptic_curves::ShortWeierstrassCurve::<F17>::new(
        F17::from_i64(5),
        F17::from_i64(7),
    )
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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn derivative_is_linear_over_same_curve_samples(
        case in arb_short_weierstrass_function_pair_case::<crate::fields::Fp17>(
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
        case in arb_short_weierstrass_function_case::<crate::fields::Fp17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        let constant = ShortWeierstrassFunction::<F17>::from_rational_function(
            case.curve.clone(),
            RationalFunction::<F17>::constant(F17::from_i64(9)),
        );

        prop_assert!(constant.derivative().is_zero());
    }

    #[test]
    fn derivative_commutes_with_conjugation(
        case in arb_short_weierstrass_function_case::<crate::fields::Fp17>(
            CurveStrategyConfig::default(),
            PolynomialStrategyConfig::default(),
        )
    ) {
        let left = case.function.derivative().conjugate();
        let right = case.function.conjugate().derivative();

        prop_assert_eq!(left, right);
    }
}
