use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::{Fp, traits::Field};
use crate::polynomials::DensePolynomial;

type F17 = Fp<17>;

#[test]
fn form_variants_report_their_shape_honestly() {
    let polynomial = DensePolynomial::<F17>::constant(F17::one());
    let in_x = DivisionPolynomialForm::x_polynomial(polynomial.clone());
    let y_times = DivisionPolynomialForm::y_times_x_polynomial(polynomial.clone());

    assert!(in_x.is_x_polynomial());
    assert!(!in_x.is_y_times_x_polynomial());
    assert_eq!(in_x.x_factor(), &polynomial);

    assert!(y_times.is_y_times_x_polynomial());
    assert!(!y_times.is_x_polynomial());
    assert_eq!(y_times.x_factor(), &polynomial);
}

#[test]
fn base_division_polynomials_cover_psi_zero_through_four() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("curve should be non-singular");

    assert_eq!(
        curve.base_division_polynomial(0).unwrap(),
        DivisionPolynomialForm::x_polynomial(DensePolynomial::<F17>::new(Vec::new()))
    );
    assert_eq!(
        curve.base_division_polynomial(1).unwrap(),
        DivisionPolynomialForm::x_polynomial(DensePolynomial::<F17>::constant(F17::one()))
    );
    assert_eq!(
        curve.base_division_polynomial(2).unwrap(),
        DivisionPolynomialForm::y_times_x_polynomial(DensePolynomial::<F17>::constant(
            F17::from_i64(2),
        ))
    );
}

#[test]
fn odd_and_even_public_polynomials_have_the_expected_shape() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("curve should be non-singular");

    assert!(matches!(
        curve.division_polynomial(3).unwrap(),
        DivisionPolynomialForm::InX(_)
    ));
    assert!(matches!(
        curve.division_polynomial(4).unwrap(),
        DivisionPolynomialForm::YTimes(_)
    ));
}

#[test]
fn zero_index_is_rejected_via_public_dispatchers() {
    let curve = ShortWeierstrassCurve::<F17>::new(F17::from_i64(2), F17::from_i64(3))
        .expect("curve should be non-singular");

    assert_eq!(
        curve.rational_x_candidates_for_division_polynomial(0),
        Err(DivisionPolynomialError::ZeroIndex)
    );
    assert_eq!(
        curve.torsion_candidates_from_division_polynomial(0),
        Err(DivisionPolynomialError::ZeroIndex)
    );
}
