use proptest::prelude::*;

use crate::elliptic_curves::analytic::{
    ApproxTolerance, ModularMatrix, UpperHalfPlanePoint,
    fundamental_domain::{FundamentalDomainReductionStatus, FundamentalDomainReductionStepReason},
};
use crate::proptest_support::{
    config::AnalyticStrategyConfig, elliptic_curves::arb_upper_half_plane_point,
};

#[test]
fn already_reduced_tau_reports_no_steps() {
    let tau = UpperHalfPlanePoint::tau_i();
    let report = tau
        .reduce_to_standard_fundamental_domain(8)
        .expect("reduction should work");

    assert_eq!(report.original_tau(), &tau);
    assert_eq!(report.reduced_tau(), &tau);
    assert_eq!(
        report.status(),
        FundamentalDomainReductionStatus::AlreadyReduced
    );
    assert!(report.steps().is_empty());
    assert_eq!(report.accumulated_matrix(), ModularMatrix::identity());
}

#[test]
fn translation_step_is_recorded_when_real_part_starts_outside_strip() {
    let tau = UpperHalfPlanePoint::from_re_im(1.7, 1.2).unwrap();
    let report = tau
        .reduce_to_standard_fundamental_domain(8)
        .expect("reduction should work");

    assert!(report.is_reduced());
    assert_eq!(
        report.steps()[0].reason(),
        FundamentalDomainReductionStepReason::RealPartOutsideCenteredStrip
    );
}

#[test]
fn inversion_step_is_recorded_when_norm_is_below_one() {
    let tau = UpperHalfPlanePoint::from_re_im(0.1, 0.8).unwrap();
    let report = tau
        .reduce_to_standard_fundamental_domain(8)
        .expect("reduction should work");

    assert!(report.is_reduced());
    assert!(
        report
            .steps()
            .iter()
            .any(|step| { step.reason() == FundamentalDomainReductionStepReason::NormLessThanOne })
    );
}

#[test]
fn step_limit_is_reported_honestly() {
    let tau = UpperHalfPlanePoint::from_re_im(3.2, 0.4).unwrap();
    let report = tau
        .reduce_to_standard_fundamental_domain(0)
        .expect("reduction should work");

    assert_eq!(
        report.status(),
        FundamentalDomainReductionStatus::StepLimitReached
    );
    assert!(!report.is_reduced());
}

#[test]
fn domain_membership_test_matches_standard_examples() {
    assert!(
        UpperHalfPlanePoint::tau_i().is_in_standard_fundamental_domain(ApproxTolerance::strict())
    );
    assert!(
        UpperHalfPlanePoint::tau_rho().is_in_standard_fundamental_domain(ApproxTolerance::strict())
    );
    assert!(
        !UpperHalfPlanePoint::from_re_im(0.8, 1.0)
            .unwrap()
            .is_in_standard_fundamental_domain(ApproxTolerance::strict())
    );
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn generic_reduction_outputs_something_in_domain_or_hits_limit(
        tau in arb_upper_half_plane_point(AnalyticStrategyConfig {
            max_real_part: 4.0,
            min_imaginary_part: 0.3,
            max_imaginary_part: 4.0,
        }),
    ) {
        let report = tau
            .reduce_to_standard_fundamental_domain(12)
            .expect("reduction should work");

        if report.is_reduced() {
            prop_assert!(report
                .reduced_tau()
                .is_in_standard_fundamental_domain(ApproxTolerance::strict()));
        } else {
            prop_assert_eq!(
                report.status(),
                FundamentalDomainReductionStatus::StepLimitReached
            );
        }
    }
}
