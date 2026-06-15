use num_complex::Complex64;
use proptest::prelude::*;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, UpperHalfPlanePoint,
    elliptic_functions::{
        EllipticFunctionTruncation, WeierstrassPDerivativeApprox,
        traits::EllipticFunctionApproximation,
    },
};
use crate::fields::ComplexApprox;
use crate::proptest_support::elliptic_curves::arb_interior_fundamental_coordinate;

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
        lattice.weierstrass_p(c(0.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
    assert_eq!(
        lattice.weierstrass_p(c(1.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
}

#[test]
fn weierstrass_p_reports_input_truncation_and_terms_used() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.3, 0.2);

    let approximation = lattice.weierstrass_p(z, truncation).unwrap();

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

    let original = lattice.weierstrass_p(z, truncation).unwrap();
    let translated = lattice.weierstrass_p(z + c(1.0, 0.0), truncation).unwrap();

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
    let approximation = lattice.weierstrass_p(c(0.5, 0.0), truncation).unwrap();

    assert!(approximation.value().im.abs() <= ComplexApprox::default_tolerance().absolute);
}

#[test]
fn weierstrass_p_derivative_rejects_lattice_points_as_poles() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();

    assert_eq!(
        lattice.weierstrass_p_derivative(c(0.0, 0.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
    assert_eq!(
        lattice.weierstrass_p_derivative(c(0.0, 1.0), truncation),
        Err(AnalyticCurveError::PointTooCloseToLatticePoint)
    );
}

#[test]
fn weierstrass_p_derivative_reports_input_truncation_and_terms_used() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.2, 0.35);

    let approximation = lattice.weierstrass_p_derivative(z, truncation).unwrap();

    assert_eq!(
        approximation,
        WeierstrassPDerivativeApprox::from_parts(
            z,
            *approximation.value(),
            truncation,
            24,
            approximation.pole_distance(),
        )
    );
}

#[test]
fn weierstrass_p_derivative_is_periodic_under_lattice_translation() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::default_educational();
    let z = c(0.18, 0.27);

    let original = lattice.weierstrass_p_derivative(z, truncation).unwrap();
    let translated = lattice
        .weierstrass_p_derivative(z + c(0.0, 1.0), truncation)
        .unwrap();

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

    let derivative = lattice.weierstrass_p_derivative(z, truncation).unwrap();
    let forward = lattice.weierstrass_p(z + c(h, 0.0), truncation).unwrap();
    let backward = lattice.weierstrass_p(z - c(h, 0.0), truncation).unwrap();
    let finite_difference = (*forward.value() - *backward.value()) / (2.0 * h);

    assert!(ComplexApprox::eq_with_tolerance(
        derivative.value(),
        &finite_difference,
        crate::numerics::ApproxTolerance::new(1.0e-5, 1.0e-5)
    ));
}

#[test]
fn weierstrass_p_is_even_approximately() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::new(14).unwrap();
    let z = c(0.23, 0.19);

    let positive = lattice.weierstrass_p(z, truncation).unwrap();
    let negative = lattice.weierstrass_p(-z, truncation).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        positive.value(),
        negative.value(),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2)
    ));
}

#[test]
fn weierstrass_p_prime_is_odd_approximately() {
    let lattice = square_lattice();
    let truncation = EllipticFunctionTruncation::new(14).unwrap();
    let z = c(0.23, 0.19);

    let positive = lattice.weierstrass_p_derivative(z, truncation).unwrap();
    let negative = lattice.weierstrass_p_derivative(-z, truncation).unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        positive.value(),
        &(-*negative.value()),
        crate::numerics::ApproxTolerance::new(1.0e-2, 1.0e-2)
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
    let p = lattice.weierstrass_p(z, truncation).unwrap();
    let dp = lattice.weierstrass_p_derivative(z, truncation).unwrap();

    assert_shared_shape(&p, z, 2, 24);
    assert_shared_shape(&dp, z, 2, 24);
    assert!(p.pole_distance() > 0.0);
    assert!(dp.pole_distance() > 0.0);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn generic_points_away_from_poles_have_finite_weierstrass_reports(
        coord in arb_interior_fundamental_coordinate(),
    ) {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();
        let z = lattice.point_from_fundamental_coordinates(coord);
        let p = lattice.weierstrass_p(z, truncation).unwrap();
        let dp = lattice.weierstrass_p_derivative(z, truncation).unwrap();

        prop_assert_eq!(*p.z(), z);
        prop_assert_eq!(*dp.z(), z);
        prop_assert_eq!(p.terms_used(), 24);
        prop_assert_eq!(dp.terms_used(), 24);
        prop_assert!(p.pole_distance() > 0.0);
        prop_assert!(dp.pole_distance() > 0.0);
        prop_assert!(p.value().re.is_finite() && p.value().im.is_finite());
        prop_assert!(dp.value().re.is_finite() && dp.value().im.is_finite());
    }

    #[test]
    fn weierstrass_functions_are_periodic_under_small_integer_lattice_shifts(
        coord in arb_interior_fundamental_coordinate(),
        m in -2i32..=2,
        n in -2i32..=2,
    ) {
        let lattice = square_lattice();
        let truncation = EllipticFunctionTruncation::default_educational();
        let tolerance = crate::numerics::ApproxTolerance::loose();
        let z = lattice.point_from_fundamental_coordinates(coord);
        let shift = c(m as f64, n as f64);

        let p = lattice.weierstrass_p(z, truncation).unwrap();
        let p_shifted = lattice.weierstrass_p(z + shift, truncation).unwrap();
        let dp = lattice.weierstrass_p_derivative(z, truncation).unwrap();
        let dp_shifted = lattice.weierstrass_p_derivative(z + shift, truncation).unwrap();

        prop_assert!(ComplexApprox::eq_with_tolerance(p.value(), p_shifted.value(), tolerance));
        prop_assert!(ComplexApprox::eq_with_tolerance(dp.value(), dp_shifted.value(), tolerance));
    }
}
