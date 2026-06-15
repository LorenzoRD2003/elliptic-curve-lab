use num_complex::Complex64;
use proptest::prelude::*;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, ComplexLattice, LatticeSumTruncation,
    UpperHalfPlanePoint,
    invariants::{ComplexAnalyticCurveLabReport, SpecialJKind, SpecialTauKind},
};
use crate::fields::{complex_approx::ComplexApprox, traits::Field};
use crate::proptest_support::{
    config::{AnalyticStrategyConfig, FieldStrategyConfig},
    elliptic_curves::arb_upper_half_plane_point,
    fields::arb_complex_approx,
};

fn standard_square_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_i())
}

fn standard_hexagonal_lattice() -> ComplexLattice {
    ComplexLattice::from_tau(UpperHalfPlanePoint::tau_rho())
}

#[test]
fn analytic_g2_matches_sixty_times_g4() {
    let lattice = standard_square_lattice();
    let truncation = LatticeSumTruncation::default_educational();

    let g2 = lattice.analytic_g2(truncation).unwrap();
    let g4 = lattice.g4_sum(truncation).unwrap();

    assert!(ComplexApprox::eq(
        &g2,
        &(Complex64::new(60.0, 0.0) * *g4.value())
    ));
}

#[test]
fn analytic_g3_matches_one_hundred_forty_times_g6() {
    let lattice = standard_hexagonal_lattice();
    let truncation = LatticeSumTruncation::default_educational();

    let g3 = lattice.analytic_g3(truncation).unwrap();
    let g6 = lattice.g6_sum(truncation).unwrap();

    assert!(ComplexApprox::eq(
        &g3,
        &(Complex64::new(140.0, 0.0) * *g6.value())
    ));
}

#[test]
fn discriminant_helper_uses_the_classical_formula() {
    let g2 = Complex64::new(3.0, -1.0);
    let g3 = Complex64::new(-2.0, 4.0);

    let discriminant = AnalyticInvariants::discriminant_from_g2_g3(&g2, &g3);
    let expected = g2.powu(3) - Complex64::new(27.0, 0.0) * g3.powu(2);

    assert_eq!(discriminant, expected);
}

#[test]
fn j_invariant_helper_rejects_nearly_singular_input() {
    let g2 = Complex64::new(0.0, 0.0);
    let g3 = Complex64::new(0.0, 0.0);

    assert_eq!(
        AnalyticInvariants::j_invariant_from_g2_g3(&g2, &g3),
        Err(AnalyticCurveError::NearlySingularAnalyticCurve)
    );
}

#[test]
fn invariants_from_tau_match_invariants_from_the_associated_lattice() {
    let tau = UpperHalfPlanePoint::tau_rho();
    let truncation = LatticeSumTruncation::larger_for_comparison();

    let from_tau = tau.analytic_invariants(truncation).unwrap();
    let lattice = ComplexLattice::from_tau(tau);
    let from_lattice = lattice.analytic_invariants(truncation).unwrap();

    assert_eq!(from_tau, from_lattice);
}

#[test]
fn lab_report_keeps_j_consistent_across_model_surfaces() {
    let report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_generic_example(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert!(ComplexApprox::eq_with_tolerance(
        report.invariants().j_invariant(),
        &report.analytic_curve().j_invariant().unwrap(),
        crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
    ));
    assert!(ComplexApprox::eq_with_tolerance(
        report.invariants().j_invariant(),
        report.short_model().j_invariant(),
        crate::numerics::ApproxTolerance::new(1.0e-6, 1.0e-6),
    ));
}

#[test]
fn special_classifications_distinguish_classical_cases() {
    let tau_i_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_i(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let tau_rho_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_rho(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();
    let generic_report = ComplexAnalyticCurveLabReport::from_tau(
        UpperHalfPlanePoint::tau_generic_example(),
        LatticeSumTruncation::new(12).unwrap(),
    )
    .unwrap();

    assert_eq!(tau_i_report.special_tau_kind(), SpecialTauKind::TauI);
    assert_eq!(tau_rho_report.special_tau_kind(), SpecialTauKind::TauRho);
    assert_eq!(generic_report.special_tau_kind(), SpecialTauKind::Generic);
    assert_eq!(tau_i_report.special_j_kind(), SpecialJKind::Near1728);
    assert_eq!(tau_rho_report.special_j_kind(), SpecialJKind::NearZero);
    assert_eq!(generic_report.special_j_kind(), SpecialJKind::Generic);
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20))]

    #[test]
    fn discriminant_formula_holds_for_generic_complex_inputs(
        g2 in arb_complex_approx(FieldStrategyConfig {
            max_real_norm: 5.0,
            max_imaginary_norm: 5.0,
            ..FieldStrategyConfig::default()
        }),
        g3 in arb_complex_approx(FieldStrategyConfig {
            max_real_norm: 5.0,
            max_imaginary_norm: 5.0,
            ..FieldStrategyConfig::default()
        }),
    ) {
        let expected = g2.powu(3) - Complex64::new(27.0, 0.0) * g3.powu(2);
        prop_assert_eq!(AnalyticInvariants::discriminant_from_g2_g3(&g2, &g3), expected);
    }

    #[test]
    fn invariants_from_tau_match_lattice_route_for_generic_tau(
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig {
            max_real_part: 0.45,
            min_imaginary_part: 0.8,
            max_imaginary_part: 2.2,
        }),
    ) {
        let truncation = LatticeSumTruncation::larger_for_comparison();
        let from_tau = tau.analytic_invariants(truncation).unwrap();
        let from_lattice = ComplexLattice::from_tau(tau).analytic_invariants(truncation).unwrap();

        prop_assert_eq!(from_tau, from_lattice);
    }
}
