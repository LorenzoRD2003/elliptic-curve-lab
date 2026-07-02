use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::traits::*;
use crate::polynomials::DensePolynomial;

type F17 = crate::fields::Fp17;
type F43 = crate::fields::Fp43;

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

#[test]
fn psi_seven_matches_the_known_sage_result_on_one_f43_curve() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::from_i64(-10), F43::from_i64(-10))
        .expect("curve should be non-singular");

    let DivisionPolynomialForm::InX(psi_seven) = curve.division_polynomial(7).unwrap() else {
        panic!("psi_7 should lie in F[x] for odd index 7");
    };

    assert_eq!(
        psi_seven.make_monic().expect("psi_7 should be non-zero"),
        DensePolynomial::new(vec![
            F43::from_i64(15),
            F43::from_i64(27),
            F43::from_i64(11),
            F43::from_i64(24),
            F43::from_i64(11),
            F43::from_i64(42),
            F43::from_i64(39),
            F43::from_i64(4),
            F43::from_i64(31),
            F43::from_i64(30),
            F43::from_i64(29),
            F43::from_i64(8),
            F43::from_i64(31),
            F43::from_i64(10),
            F43::one(),
            F43::from_i64(40),
            F43::from_i64(37),
            F43::from_i64(38),
            F43::from_i64(14),
            F43::from_i64(34),
            F43::from_i64(26),
            F43::from_i64(11),
            F43::from_i64(33),
            F43::zero(),
            F43::one(),
        ])
    );
}
