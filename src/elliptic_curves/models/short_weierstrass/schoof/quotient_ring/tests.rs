use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    short_weierstrass::schoof::{
        QuotientInverseResult, ReducedCurveFunction, ReducedCurveQuotient,
    },
};
use crate::fields::{Fp, traits::Field};
use crate::polynomials::DensePolynomial;

type F7 = Fp<7>;

fn sample_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
        .expect("sample curve should be smooth")
}

#[test]
fn quotient_rejects_the_zero_polynomial_modulus() {
    let curve = sample_curve();

    assert!(matches!(
        ReducedCurveQuotient::new(curve, DensePolynomial::new(Vec::new())),
        Err(CurveError::ZeroReducedCurveQuotientModulus)
    ));
}

#[test]
fn quotient_normalizes_the_modulus_to_monic_form() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(2), F7::from_i64(4)]),
    )
    .expect("non-zero modulus should define a quotient");

    assert_eq!(
        quotient.modulus(),
        &DensePolynomial::new(vec![F7::from_i64(4), F7::one()])
    );
    assert_eq!(quotient.curve().a(), &F7::from_i64(2));
    assert_eq!(quotient.curve().b(), &F7::from_i64(3));
}

#[test]
fn reduce_poly_returns_the_unique_remainder_modulo_g() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let polynomial = DensePolynomial::new(vec![
        F7::from_i64(3),
        F7::from_i64(1),
        F7::from_i64(5),
        F7::from_i64(2),
    ]);

    assert_eq!(
        quotient.reduce_poly(&polynomial),
        polynomial.rem(quotient.modulus()).unwrap()
    );
}

#[test]
fn distinguished_classes_use_canonical_reduced_coordinates() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");

    let zero = quotient.zero();
    let one = quotient.one();
    let x = quotient.x();
    let y = quotient.y();

    assert!(zero.is_zero());
    assert_eq!(one.x_part(), &DensePolynomial::constant(F7::one()));
    assert!(one.y_part().is_zero());
    assert_eq!(
        x.x_part(),
        &DensePolynomial::new(vec![F7::zero(), F7::one()])
    );
    assert!(x.y_part().is_zero());
    assert!(y.x_part().is_zero());
    assert_eq!(y.y_part(), &DensePolynomial::constant(F7::one()));
}

#[test]
fn reduced_curve_function_addition_and_subtraction_are_componentwise_mod_g() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let left = ReducedCurveFunction::new(
        &quotient,
        DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(1), F7::from_i64(2)]),
        DensePolynomial::new(vec![F7::from_i64(6), F7::from_i64(5)]),
    );
    let right = ReducedCurveFunction::new(
        &quotient,
        DensePolynomial::new(vec![F7::from_i64(4), F7::from_i64(2)]),
        DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(6), F7::from_i64(3)]),
    );

    let sum = left.add(&quotient, &right);
    let difference = sum.sub(&quotient, &right);

    assert_eq!(
        sum.x_part(),
        &quotient.reduce_poly(&left.x_part().add(right.x_part()))
    );
    assert_eq!(
        sum.y_part(),
        &quotient.reduce_poly(&left.y_part().add(right.y_part()))
    );
    assert_eq!(difference.x_part(), left.x_part());
    assert_eq!(difference.y_part(), left.y_part());
}

#[test]
fn reduced_curve_function_negation_produces_additive_inverse() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let value = ReducedCurveFunction::new(
        &quotient,
        DensePolynomial::new(vec![F7::from_i64(2), F7::from_i64(4)]),
        DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(1)]),
    );

    let negated = value.neg(&quotient);
    let sum = value.add(&quotient, &negated);

    assert!(sum.is_zero());
}

#[test]
fn reduced_curve_function_multiplication_uses_y_squared_equals_f_of_x() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let left = ReducedCurveFunction::new(
        &quotient,
        DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(2)]),
        DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(1)]),
    );
    let right = ReducedCurveFunction::new(
        &quotient,
        DensePolynomial::new(vec![F7::from_i64(6), F7::from_i64(5)]),
        DensePolynomial::new(vec![F7::from_i64(4)]),
    );

    let product = left.mul(&quotient, &right);
    let cubic = quotient.curve().to_cubic();
    let expected_x = quotient.reduce_poly(
        &left
            .x_part()
            .mul(right.x_part())
            .add(&cubic.mul(&left.y_part().mul(right.y_part()))),
    );
    let expected_y = quotient.reduce_poly(
        &left
            .x_part()
            .mul(right.y_part())
            .add(&right.x_part().mul(left.y_part())),
    );

    assert_eq!(product.x_part(), &expected_x);
    assert_eq!(product.y_part(), &expected_y);
}

#[test]
fn try_invert_poly_returns_a_reduced_inverse_for_units() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let unit = DensePolynomial::new(vec![F7::one(), F7::one()]);

    let QuotientInverseResult::Inverse(inverse) = quotient.try_invert_poly(&unit) else {
        panic!("x + 1 should be invertible modulo x^2 + 1 over F7");
    };

    let product = quotient.reduce_poly(&unit.mul(&inverse));
    assert_eq!(product, DensePolynomial::constant(F7::one()));
}

#[test]
fn try_invert_poly_returns_gcd_witness_for_non_units() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::from_i64(-1), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");
    let non_unit = DensePolynomial::new(vec![F7::from_i64(-1), F7::one()]);

    let QuotientInverseResult::NonUnit { witness_gcd } = quotient.try_invert_poly(&non_unit) else {
        panic!("x - 1 should not be invertible modulo x^2 - 1");
    };

    assert_eq!(
        witness_gcd,
        DensePolynomial::new(vec![F7::from_i64(-1), F7::one()])
    );
}

#[test]
fn try_invert_poly_reports_zero_as_a_non_unit_with_full_modulus_witness() {
    let curve = sample_curve();
    let quotient = ReducedCurveQuotient::new(
        curve,
        DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
    )
    .expect("non-zero modulus should define a quotient");

    let QuotientInverseResult::NonUnit { witness_gcd } =
        quotient.try_invert_poly(&DensePolynomial::new(Vec::new()))
    else {
        panic!("the zero class should never be invertible");
    };

    assert_eq!(witness_gcd, quotient.modulus().clone());
}
