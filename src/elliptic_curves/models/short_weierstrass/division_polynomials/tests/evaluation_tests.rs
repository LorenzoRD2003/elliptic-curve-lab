use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::DivisionPolynomialError, traits::AffineCurveModel,
};
use crate::fields::{Fp, traits::Field};

type F23 = Fp<23>;

#[test]
fn x_criterion_evaluation_matches_specialized_evaluators() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");
    let x = F23::from_i64(1);

    assert_eq!(
        curve
            .evaluate_division_polynomial_x_criterion(5, &x)
            .unwrap(),
        curve.evaluate_odd_division_polynomial_at_x(5, &x).unwrap()
    );
    assert_eq!(
        curve
            .evaluate_division_polynomial_x_criterion(6, &x)
            .unwrap(),
        curve
            .evaluate_even_division_polynomial_factor_at_x(6, &x)
            .unwrap()
    );
}

#[test]
fn point_evaluation_rejects_identity_and_off_curve_inputs() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");

    assert_eq!(
        curve.evaluate_division_polynomial_at_point(3, &AffinePoint::<F23>::infinity()),
        Err(DivisionPolynomialError::PointAtInfinityNotSupported)
    );
    assert_eq!(
        curve.evaluate_division_polynomial_at_point(
            3,
            &AffinePoint::<F23>::new(F23::zero(), F23::zero())
        ),
        Err(DivisionPolynomialError::Curve(
            crate::elliptic_curves::CurveError::PointNotOnCurve
        ))
    );
}

#[test]
fn odd_point_evaluation_depends_only_on_x_and_even_changes_sign() {
    let curve = ShortWeierstrassCurve::<F23>::new(F23::from_i64(2), F23::from_i64(3))
        .expect("curve should be non-singular");
    let point = curve.point(F23::from_i64(1), F23::from_i64(12)).unwrap();

    let odd = curve
        .evaluate_division_polynomial_at_point(5, &point)
        .unwrap();
    let odd_neg = curve
        .evaluate_division_polynomial_at_point(5, &point.neg())
        .unwrap();
    let even = curve
        .evaluate_division_polynomial_at_point(6, &point)
        .unwrap();
    let even_neg = curve
        .evaluate_division_polynomial_at_point(6, &point.neg())
        .unwrap();

    assert!(F23::eq(&odd, &odd_neg));
    assert!(F23::eq(&even_neg, &F23::neg(&even)));
}
