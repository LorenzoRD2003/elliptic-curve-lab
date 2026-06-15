use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, EllipticFunctionTruncation, UpperHalfPlanePoint,
    zeta::WeierstrassZetaApprox,
};
use crate::fields::ComplexApprox;

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

#[test]
fn weierstrass_zeta_rejects_lattice_points_as_poles() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();

    assert_eq!(
        lattice.weierstrass_zeta(c(0.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
    assert_eq!(
        lattice.weierstrass_zeta(c(1.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
}

#[test]
fn weierstrass_zeta_reports_input_truncation_and_terms_used() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.3, 0.2);

    let approximation = lattice.weierstrass_zeta(z, truncation).unwrap();

    assert_eq!(
        approximation,
        WeierstrassZetaApprox::new(
            z,
            *approximation.value(),
            truncation,
            24,
            approximation.pole_distance(),
        )
    );
}

#[test]
fn weierstrass_zeta_is_odd() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.31, 0.22);

    let positive = lattice.weierstrass_zeta(z, truncation).unwrap();
    let negative = lattice.weierstrass_zeta(-z, truncation).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        negative.value(),
        &(-*positive.value()),
        ComplexApprox::default_tolerance()
    ));
}

#[test]
fn weierstrass_zeta_derivative_matches_minus_weierstrass_p_by_finite_difference() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.23, 0.19);
    let h = 1.0e-6;

    let p = lattice.weierstrass_p(z, truncation).unwrap();
    let zeta_forward = lattice.weierstrass_zeta(z + c(h, 0.0), truncation).unwrap();
    let zeta_backward = lattice.weierstrass_zeta(z - c(h, 0.0), truncation).unwrap();
    let finite_difference = (*zeta_forward.value() - *zeta_backward.value()) / (2.0 * h);

    assert!(ComplexApprox::eq_with_tolerance(
        &finite_difference,
        &(-*p.value()),
        crate::numerics::ApproxTolerance::new(1.0e-5, 1.0e-5)
    ));
}
