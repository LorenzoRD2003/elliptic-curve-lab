use num_complex::Complex64;

use super::{
    EllipticFunctionApproximation, EllipticFunctionTruncation, HasPoleDistance,
    WeierstrassPDerivativeApprox, weierstrass_p, weierstrass_p_derivative,
};
use crate::{
    elliptic_curves::analytic::{AnalyticCurveError, ComplexLattice, UpperHalfPlanePoint},
    fields::ComplexApprox,
};

fn c(re: f64, im: f64) -> Complex64 {
    Complex64::new(re, im)
}

fn square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

#[test]
fn elliptic_function_truncation_requires_positive_radius() {
    assert_eq!(
        EllipticFunctionTruncation::new(0),
        Err(AnalyticCurveError::InvalidTruncationRadius)
    );
}

#[test]
fn elliptic_function_truncation_exposes_radius() {
    let truncation = EllipticFunctionTruncation::new(3).unwrap();

    assert_eq!(truncation.radius(), 3);
}

#[test]
fn educational_default_radius_is_explicit() {
    let truncation = EllipticFunctionTruncation::default_educational();

    assert_eq!(truncation.radius(), 2);
}

#[test]
fn weierstrass_p_rejects_lattice_points_as_poles() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();

    assert_eq!(
        weierstrass_p(&lattice, c(0.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
    assert_eq!(
        weierstrass_p(&lattice, c(1.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
}

#[test]
fn weierstrass_p_reports_input_truncation_and_terms_used() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.3, 0.2);

    let approximation = weierstrass_p(&lattice, z, truncation).unwrap();

    assert_eq!(*approximation.z(), z);
    assert_eq!(approximation.truncation(), truncation);
    assert_eq!(approximation.terms_used(), 24);
    assert!(approximation.pole_distance() > 0.0);
    assert!(approximation.value().re.is_finite());
    assert!(approximation.value().im.is_finite());
}

#[test]
fn weierstrass_p_is_periodic_under_lattice_translation() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.31, 0.22);

    let original = weierstrass_p(&lattice, z, truncation).unwrap();
    let translated = weierstrass_p(&lattice, z + c(1.0, 0.0), truncation).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        original.value(),
        translated.value(),
        ComplexApprox::default_tolerance()
    ));
}

#[test]
fn weierstrass_p_at_half_period_is_real_for_the_square_lattice() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let approximation = weierstrass_p(&lattice, c(0.5, 0.0), truncation).unwrap();

    assert!(approximation.value().im.abs() <= ComplexApprox::default_tolerance().absolute);
}

#[test]
fn weierstrass_p_derivative_rejects_lattice_points_as_poles() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();

    assert_eq!(
        weierstrass_p_derivative(&lattice, c(0.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
    assert_eq!(
        weierstrass_p_derivative(&lattice, c(0.0, 1.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
}

#[test]
fn weierstrass_p_derivative_reports_input_truncation_and_terms_used() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.2, 0.35);

    let approximation = weierstrass_p_derivative(&lattice, z, truncation).unwrap();

    assert_eq!(
        approximation,
        WeierstrassPDerivativeApprox {
            z,
            value: *approximation.value(),
            truncation,
            terms_used: 24,
            pole_distance: approximation.pole_distance(),
        }
    );
}

#[test]
fn weierstrass_p_derivative_is_periodic_under_lattice_translation() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.18, 0.27);

    let original = weierstrass_p_derivative(&lattice, z, truncation).unwrap();
    let translated = weierstrass_p_derivative(&lattice, z + c(0.0, 1.0), truncation).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        original.value(),
        translated.value(),
        ComplexApprox::default_tolerance()
    ));
}

#[test]
fn weierstrass_p_derivative_matches_a_centered_finite_difference_of_weierstrass_p() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.23, 0.19);
    let h = 1.0e-6;

    let derivative = weierstrass_p_derivative(&lattice, z, truncation).unwrap();
    let forward = weierstrass_p(&lattice, z + c(h, 0.0), truncation).unwrap();
    let backward = weierstrass_p(&lattice, z - c(h, 0.0), truncation).unwrap();
    let finite_difference = (*forward.value() - *backward.value()) / (2.0 * h);

    assert!(ComplexApprox::eq_with_tolerance(
        derivative.value(),
        &finite_difference,
        crate::numerics::ApproxTolerance::new(1.0e-5, 1.0e-5)
    ));
}

#[test]
fn approximation_trait_exposes_shared_report_metadata() {
    fn assert_shared_shape<T: EllipticFunctionApproximation>(
        approximation: &T,
        expected_z: Complex64,
        expected_radius: usize,
        expected_terms: usize,
    ) {
        assert_eq!(*approximation.z(), expected_z);
        assert_eq!(approximation.truncation().radius(), expected_radius);
        assert_eq!(approximation.terms_used(), expected_terms);
        assert!(approximation.value().re.is_finite());
        assert!(approximation.value().im.is_finite());
    }

    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.17, 0.29);
    let p = weierstrass_p(&lattice, z, truncation).unwrap();
    let dp = weierstrass_p_derivative(&lattice, z, truncation).unwrap();

    assert_shared_shape(&p, z, 2, 24);
    assert_shared_shape(&dp, z, 2, 24);
    assert!(p.pole_distance() > 0.0);
    assert!(dp.pole_distance() > 0.0);
}
