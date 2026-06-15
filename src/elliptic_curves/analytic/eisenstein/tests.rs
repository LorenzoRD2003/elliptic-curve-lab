use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, LatticeSumTruncation, UpperHalfPlanePoint,
    eisenstein::{EisensteinSumApprox, TruncationConvergenceReport},
};
use crate::fields::complex_approx::ComplexApprox;

fn standard_square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

fn standard_hexagonal_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_rho())
}

#[test]
fn eisenstein_sum_rejects_non_convergent_low_weights() {
    let lattice = standard_square_lattice();
    let truncation = LatticeSumTruncation::default_educational();

    assert_eq!(
        lattice.eisenstein_sum(0, truncation),
        Err(AnalyticCurveError::InvalidEisensteinWeight)
    );
    assert_eq!(
        lattice.eisenstein_sum(1, truncation),
        Err(AnalyticCurveError::InvalidEisensteinWeight)
    );
    assert_eq!(
        lattice.eisenstein_sum(2, truncation),
        Err(AnalyticCurveError::InvalidEisensteinWeight)
    );
}

#[test]
fn odd_weight_eisenstein_sum_short_circuits_to_exact_zero() {
    let lattice = standard_hexagonal_lattice();
    let truncation = LatticeSumTruncation::larger_for_comparison();

    let sum = lattice.eisenstein_sum(5, truncation).unwrap();

    assert_eq!(
        sum,
        EisensteinSumApprox::new(
            5,
            Complex64::new(0.0, 0.0),
            truncation,
            truncation.nonzero_terms_in_square_box()
        )
    );
}

#[test]
fn g4_sum_matches_direct_even_weight_dispatch() {
    let lattice = standard_square_lattice();
    let truncation = LatticeSumTruncation::default_educational();

    let generic = lattice.eisenstein_sum(4, truncation).unwrap();
    let specialized = lattice.g4_sum(truncation).unwrap();

    assert_eq!(generic, specialized);
}

#[test]
fn g6_sum_matches_direct_even_weight_dispatch() {
    let lattice = standard_square_lattice();
    let truncation = LatticeSumTruncation::default_educational();

    let generic = lattice.eisenstein_sum(6, truncation).unwrap();
    let specialized = lattice.g6_sum(truncation).unwrap();

    assert_eq!(generic, specialized);
}

#[test]
fn larger_truncation_comparison_rejects_non_increasing_pairs() {
    let lattice = standard_square_lattice();
    let same = LatticeSumTruncation::default_educational();
    let larger = LatticeSumTruncation::larger_for_comparison();

    assert_eq!(
        lattice.compare_eisenstein_truncations(4, same, same),
        Err(AnalyticCurveError::InvalidTruncationComparison)
    );
    assert_eq!(
        lattice.compare_eisenstein_truncations(4, larger, same),
        Err(AnalyticCurveError::InvalidTruncationComparison)
    );
}

#[test]
fn truncation_comparison_records_both_values_and_difference() {
    let lattice = standard_square_lattice();
    let small = LatticeSumTruncation::default_educational();
    let large = LatticeSumTruncation::larger_for_comparison();

    let report = lattice
        .compare_eisenstein_truncations(4, small, large)
        .unwrap();
    let small_sum = lattice.g4_sum(small).unwrap();
    let large_sum = lattice.g4_sum(large).unwrap();

    assert_eq!(
        report,
        TruncationConvergenceReport::new(
            small,
            large,
            crate::numerics::ComplexDifferenceReport::new(*large_sum.value(), *small_sum.value())
        )
    );
    assert!(report.absolute_difference().is_finite());
}

#[test]
fn larger_boxes_refine_square_lattice_weight_four_sum() {
    let lattice = standard_square_lattice();
    let small = lattice
        .g4_sum(LatticeSumTruncation::default_educational())
        .unwrap();
    let larger = lattice
        .g4_sum(LatticeSumTruncation::larger_for_comparison())
        .unwrap();

    assert!(small.value().re.is_finite());
    assert!(larger.value().re.is_finite());
    assert!(ComplexApprox::eq_with_tolerance(
        small.value(),
        larger.value(),
        crate::numerics::ApproxTolerance::new(1.0e-1, 1.0e-1),
    ));
}
