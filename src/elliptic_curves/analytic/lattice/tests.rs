use num_complex::Complex64;

use super::{
    ComplexLattice, ComplexTorusPoint, FundamentalParallelogramCoordinate, LatticeIndexPoint,
};
use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};
use crate::fields::ComplexApprox;

#[test]
fn constructor_accepts_positive_oriented_basis() {
    let lattice = ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0))
        .expect("standard basis should be valid");

    assert_eq!(lattice.omega1(), &Complex64::new(1.0, 0.0));
    assert_eq!(lattice.omega2(), &Complex64::new(0.0, 1.0));
}

#[test]
fn constructor_rejects_degenerate_basis() {
    assert_eq!(
        ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(2.0, 0.0)),
        Err(AnalyticCurveError::DegenerateLattice)
    );
}

#[test]
fn constructor_rejects_negative_orientation() {
    assert_eq!(
        ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(0.0, -1.0)),
        Err(AnalyticCurveError::NonPositiveLatticeOrientation)
    );
}

#[test]
fn from_tau_builds_standard_z_plus_z_tau_lattice() {
    let tau = UpperHalfPlanePoint::tau_i();
    let lattice = ComplexLattice::from_tau(tau);

    assert_eq!(lattice.omega1(), &Complex64::new(1.0, 0.0));
    assert_eq!(lattice.omega2(), &Complex64::new(0.0, 1.0));
}

#[test]
fn tau_recovers_the_expected_ratio() {
    let lattice = ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(-0.5, 2.0))
        .expect("basis should be valid");

    let tau = lattice.tau().expect("ratio should lie in upper half-plane");

    assert_eq!(tau.tau(), &Complex64::new(-0.5, 2.0));
}

#[test]
fn lattice_point_forms_integer_linear_combinations() {
    let lattice = ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(0.0, 2.0))
        .expect("basis should be valid");

    assert_eq!(lattice.lattice_point(3, -1), Complex64::new(3.0, -2.0));
}

#[test]
fn oriented_area_and_covolume_match_for_positive_basis() {
    let lattice = ComplexLattice::new(Complex64::new(2.0, 0.0), Complex64::new(1.0, 3.0))
        .expect("basis should be valid");

    assert_eq!(lattice.oriented_area(), 6.0);
    assert_eq!(lattice.covolume(), 6.0);
}

#[test]
fn validate_basis_matches_constructor_rules() {
    assert_eq!(
        ComplexLattice::validate_basis(Complex64::new(1.0, 0.0), Complex64::new(0.5, 1.0)),
        Ok(())
    );
    assert_eq!(
        ComplexLattice::validate_basis(Complex64::new(1.0, 0.0), Complex64::new(0.5, -1.0)),
        Err(AnalyticCurveError::NonPositiveLatticeOrientation)
    );
}

#[test]
fn validate_basis_distinguishes_nearly_degenerate_from_negative_orientation() {
    let tolerance = ComplexApprox::default_tolerance().absolute;

    assert_eq!(
        ComplexLattice::validate_basis(
            Complex64::new(1.0, 0.0),
            Complex64::new(0.5, 0.25 * tolerance)
        ),
        Err(AnalyticCurveError::DegenerateLattice)
    );
    assert_eq!(
        ComplexLattice::validate_basis(Complex64::new(1.0, 0.0), Complex64::new(0.5, -1.0)),
        Err(AnalyticCurveError::NonPositiveLatticeOrientation)
    );
}

#[test]
fn lattice_points_in_box_includes_origin_and_uses_lexicographic_order() {
    let lattice = ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0))
        .expect("basis should be valid");

    let points = lattice.lattice_points_in_box(1);

    assert_eq!(points.len(), 9);
    assert_eq!(
        points.first(),
        Some(&LatticeIndexPoint {
            m: -1,
            n: -1,
            value: Complex64::new(-1.0, -1.0),
        })
    );
    assert_eq!(
        points[4],
        LatticeIndexPoint {
            m: 0,
            n: 0,
            value: Complex64::new(0.0, 0.0),
        }
    );
    assert_eq!(
        points.last(),
        Some(&LatticeIndexPoint {
            m: 1,
            n: 1,
            value: Complex64::new(1.0, 1.0),
        })
    );
}

#[test]
fn nonzero_lattice_points_in_box_omits_only_the_origin() {
    let lattice = ComplexLattice::new(Complex64::new(1.0, 0.0), Complex64::new(0.0, 2.0))
        .expect("basis should be valid");

    let points = lattice.nonzero_lattice_points_in_box(1);

    assert_eq!(points.len(), 8);
    assert!(!points.iter().any(|point| point.m == 0 && point.n == 0));
    assert!(points.iter().any(|point| {
        point
            == &LatticeIndexPoint {
                m: 0,
                n: 1,
                value: Complex64::new(0.0, 2.0),
            }
    }));
}

#[test]
fn coordinate_constructor_accepts_exactly_the_half_open_unit_square() {
    let coordinate = FundamentalParallelogramCoordinate::new(0.0, 0.75)
        .expect("coordinates inside the canonical region should be valid");

    assert!(coordinate.is_in_half_open_unit_square());
    assert_eq!(coordinate.u(), 0.0);
    assert_eq!(coordinate.v(), 0.75);
    assert_eq!(
        FundamentalParallelogramCoordinate::new(1.0, 0.75),
        Err(AnalyticCurveError::PointNotInFundamentalParallelogram)
    );
    assert_eq!(
        FundamentalParallelogramCoordinate::new(-0.1, 0.75),
        Err(AnalyticCurveError::PointNotInFundamentalParallelogram)
    );
    assert_eq!(
        FundamentalParallelogramCoordinate::new(f64::NAN, 0.75),
        Err(AnalyticCurveError::NumericalComparisonFailed)
    );
}

#[test]
fn coordinate_reduction_snaps_near_boundary_values_back_to_zero() {
    let tolerance = ComplexApprox::default_tolerance();
    let coordinate = FundamentalParallelogramCoordinate::reduce_mod_unit_square(
        -0.25 * tolerance.absolute,
        1.0 + 0.25 * tolerance.absolute,
    )
    .expect("reduction should land inside the canonical region");

    assert_eq!(
        coordinate,
        FundamentalParallelogramCoordinate::new(0.0, 0.0)
            .expect("the snapped origin should be valid")
    );
}

#[test]
fn point_from_fundamental_coordinates_forms_the_expected_representative() {
    let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
    let point = lattice.point_from_fundamental_coordinates(
        FundamentalParallelogramCoordinate::new(0.25, 0.5)
            .expect("canonical coordinates should be valid"),
    );

    assert_eq!(point, Complex64::new(0.25, 0.5));
}

#[test]
fn reducing_a_complex_point_recovers_canonical_standard_lattice_coordinates() {
    let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
    let coordinate = lattice
        .reduce_complex_point_to_fundamental_coordinates(Complex64::new(2.25, -1.25))
        .expect("standard reduction should succeed");

    assert_eq!(
        coordinate,
        FundamentalParallelogramCoordinate::new(0.25, 0.75)
            .expect("canonical coordinates should be valid")
    );
}

#[test]
fn reducing_a_complex_point_roundtrips_for_a_generic_lattice() {
    let lattice = ComplexLattice::new(Complex64::new(2.0, 0.0), Complex64::new(1.0, 3.0))
        .expect("basis should be valid");
    let representative = lattice.point_from_fundamental_coordinates(
        FundamentalParallelogramCoordinate::reduce_mod_unit_square(1.25, -0.1)
            .expect("reduction should produce a canonical coordinate"),
    );
    let reduced = lattice
        .reduce_complex_point_to_fundamental_coordinates(representative)
        .expect("generic reduction should succeed");
    let tolerance = ComplexApprox::default_tolerance();

    assert!(tolerance.real_close(reduced.u(), 0.25));
    assert!(tolerance.real_close(reduced.v(), 0.9));
    assert!(reduced.is_in_half_open_unit_square());
}

#[test]
fn torus_point_stores_the_canonical_reduced_coordinate() {
    let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
    let torus_point = lattice
        .reduce_complex_point_to_torus_point(Complex64::new(-0.25, 1.5))
        .expect("torus reduction should succeed");

    let expected = ComplexTorusPoint::new(
        FundamentalParallelogramCoordinate::new(0.75, 0.5)
            .expect("canonical coordinates should be valid"),
    );

    assert!(lattice.torus_points_eq(&torus_point, &expected));
    assert_eq!(
        torus_point.coordinate(),
        &FundamentalParallelogramCoordinate::new(0.75, 0.5)
            .expect("canonical coordinates should be valid")
    );
}

#[test]
fn torus_points_eq_uses_explicit_lattice_relative_equality() {
    let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());
    let lhs = lattice
        .reduce_complex_point_to_torus_point(Complex64::new(0.25, 1.5))
        .expect("torus reduction should succeed");
    let rhs = lattice
        .reduce_complex_point_to_torus_point(Complex64::new(1.25, 0.5))
        .expect("torus reduction should succeed");
    let different = lattice
        .reduce_complex_point_to_torus_point(Complex64::new(0.25, 0.75))
        .expect("torus reduction should succeed");

    assert!(lattice.torus_points_eq(&lhs, &rhs));
    assert!(!lattice.torus_points_eq(&lhs, &different));
}

#[test]
fn reducing_non_finite_complex_points_fails_honestly() {
    let lattice = ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i());

    assert_eq!(
        lattice.reduce_complex_point_to_fundamental_coordinates(Complex64::new(f64::NAN, 0.0)),
        Err(AnalyticCurveError::NumericalComparisonFailed)
    );
    assert!(matches!(
        lattice.reduce_complex_point_to_torus_point(Complex64::new(f64::INFINITY, 0.0)),
        Err(AnalyticCurveError::NumericalComparisonFailed)
    ));
}
