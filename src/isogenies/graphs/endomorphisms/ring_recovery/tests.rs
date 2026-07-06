use num_bigint::BigUint;

use crate::elliptic_curves::{
    ShortWeierstrassCurve, endomorphisms::candidate_sets::EndomorphismRingCandidateSet,
    endomorphisms::quadratic_orders::QuadraticDiscriminant,
};
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId,
    endomorphisms::{
        EndomorphismRingLevelRecoveryError, EndomorphismRingLevelRecoveryReport,
        LocalEndomorphismRingLevelReport, ShortestFloorPathReport,
    },
};

type F41 = crate::fields::Fp41;
type F17 = crate::fields::Fp17;
type Curve41 = ShortWeierstrassCurve<F41>;
type Curve17 = ShortWeierstrassCurve<F17>;

fn candidate_set(discriminant: i64) -> EndomorphismRingCandidateSet {
    EndomorphismRingCandidateSet::from_discriminant(&QuadraticDiscriminant::new(discriminant))
        .expect("test discriminant should produce imaginary quadratic candidates")
}

fn local_report_from_discriminant_and_path(
    discriminant: i64,
    path: Vec<IsogenyGraphNodeId>,
) -> Result<LocalEndomorphismRingLevelReport, EndomorphismRingLevelRecoveryError> {
    local_report_at_from_discriminant_and_path(BigUint::from(2u8), discriminant, path)
}

fn local_report_at_from_discriminant_and_path(
    prime: BigUint,
    discriminant: i64,
    path: Vec<IsogenyGraphNodeId>,
) -> Result<LocalEndomorphismRingLevelReport, EndomorphismRingLevelRecoveryError> {
    let local_view = candidate_set(discriminant)
        .local_view_at(&prime)
        .expect("test prime should be a valid local prime");
    let start = path.first().copied().expect("test path should be nonempty");
    let floor = path.last().copied().expect("test path should be nonempty");
    let floor_path = ShortestFloorPathReport::new(prime, start, floor, path);

    LocalEndomorphismRingLevelReport::from_local_view_and_floor_path(local_view, floor_path)
}

fn f41_floor_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn f17_root_recovery_curve() -> Curve17 {
    Curve17::new(F17::from_i64(11), F17::from_i64(5)).expect("valid F17 curve")
}

#[test]
fn local_report_recovers_conductor_valuation_by_subtracting_floor_distance() {
    let report = local_report_from_discriminant_and_path(
        -64,
        vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)],
    )
    .expect("δ <= e should produce a local recovery report");

    assert_eq!(report.node_id(), IsogenyGraphNodeId(0));
    assert_eq!(report.prime(), &BigUint::from(2u8));
    assert_eq!(report.frobenius_conductor_valuation(), 2);
    assert_eq!(report.distance_to_floor(), 1);
    assert_eq!(report.recovered_conductor_valuation(), 1);
    assert_eq!(
        report.floor_path().path(),
        &[IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]
    );
}

#[test]
fn local_report_rejects_floor_distance_larger_than_frobenius_conductor_level() {
    let error = local_report_from_discriminant_and_path(
        -16,
        vec![
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2),
        ],
    )
    .expect_err("δ > e should be rejected as incompatible local evidence");

    assert_eq!(
        error,
        EndomorphismRingLevelRecoveryError::DistanceExceedsFrobeniusConductorValuation {
            node_id: IsogenyGraphNodeId(0),
            prime: BigUint::from(2u8),
            distance_to_floor: 2,
            frobenius_conductor_valuation: 1,
        }
    );
}

#[test]
fn graph_recovers_trivial_local_level_when_start_is_already_on_floor() {
    let graph = IsogenyGraphBuilder::new(f41_floor_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");

    let report = graph
        .recover_endomorphism_ring_level_at(IsogenyGraphNodeId(0), &BigUint::from(5u8))
        .expect("floor node should recover the local exponent with δ = 0");

    assert_eq!(report.node_id(), IsogenyGraphNodeId(0));
    assert_eq!(report.prime(), &BigUint::from(5u8));
    assert_eq!(report.distance_to_floor(), 0);
    assert_eq!(
        report.recovered_conductor_valuation(),
        report.frobenius_conductor_valuation()
    );
}

#[test]
fn graph_recovers_global_ring_report_from_supplied_local_primes() {
    let graph = IsogenyGraphBuilder::new(f17_root_recovery_curve(), 2)
        .max_depth(0)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("root graph should build");

    let report = graph
        .recover_endomorphism_ring_at(IsogenyGraphNodeId(0), &[BigUint::from(2u8)])
        .expect("official graph recovery should build local evidence and assemble it");

    assert_eq!(report.node_id(), Some(IsogenyGraphNodeId(0)));
    assert!(report.is_complete());
    assert_eq!(report.local_reports().len(), 1);
    assert_eq!(report.local_reports()[0].node_id(), IsogenyGraphNodeId(0));
    assert_eq!(report.local_reports()[0].prime(), &BigUint::from(2u8));
    assert!(report.recovered_conductor().is_some());
    assert!(report.recovered_order().is_some());
}

#[test]
fn global_report_recovers_order_when_all_prime_factors_are_covered() {
    let candidate_set = candidate_set(-108);
    let two_report = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("2-local report should build");
    let three_report = local_report_at_from_discriminant_and_path(
        BigUint::from(3u8),
        -108,
        vec![IsogenyGraphNodeId(0), IsogenyGraphNodeId(2)],
    )
    .expect("3-local report should build");

    let report = EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set.clone(),
        vec![three_report, two_report],
    )
    .expect("complete local evidence should assemble");

    assert_eq!(report.node_id(), Some(IsogenyGraphNodeId(0)));
    assert_eq!(report.candidate_set(), &candidate_set);
    assert!(report.is_complete());
    assert!(report.missing_primes().is_empty());
    assert_eq!(report.local_reports()[0].prime(), &BigUint::from(2u8));
    assert_eq!(report.local_reports()[1].prime(), &BigUint::from(3u8));
    assert_eq!(report.recovered_conductor(), Some(&BigUint::from(2u8)));
    assert_eq!(
        report
            .recovered_order()
            .expect("complete evidence should recover an order")
            .conductor(),
        &BigUint::from(2u8)
    );
}

#[test]
fn global_report_stays_partial_when_a_prime_factor_is_missing() {
    let candidate_set = candidate_set(-108);
    let two_report = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("2-local report should build");

    let report =
        EndomorphismRingLevelRecoveryReport::from_local_reports(candidate_set, vec![two_report])
            .expect("partial local evidence should still assemble");

    assert!(!report.is_complete());
    assert_eq!(report.missing_primes(), &[BigUint::from(3u8)]);
    assert!(report.recovered_conductor().is_none());
    assert!(report.recovered_order().is_none());
}

#[test]
fn global_report_recovers_maximal_order_when_frobenius_conductor_is_one() {
    let candidate_set = candidate_set(-3);

    let report = EndomorphismRingLevelRecoveryReport::from_local_reports(candidate_set, vec![])
        .expect("no local reports are needed when v = 1");

    assert_eq!(report.node_id(), None);
    assert!(report.is_complete());
    assert!(report.missing_primes().is_empty());
    assert_eq!(report.recovered_conductor(), Some(&BigUint::from(1u8)));
    assert!(
        report
            .recovered_order()
            .expect("v = 1 should recover O_K")
            .is_maximal()
    );
}

#[test]
fn global_report_rejects_duplicate_local_primes() {
    let first = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("first 2-local report should build");
    let second = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("second 2-local report should build");

    let error = EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set(-108),
        vec![first, second],
    )
    .expect_err("duplicated local primes should be rejected");

    assert_eq!(
        error,
        EndomorphismRingLevelRecoveryError::DuplicateLocalPrime {
            prime: BigUint::from(2u8),
        }
    );
}

#[test]
fn global_report_rejects_prime_not_dividing_frobenius_conductor() {
    let report_for_three = local_report_at_from_discriminant_and_path(
        BigUint::from(3u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("3-local report should build");

    let error = EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set(-16),
        vec![report_for_three],
    )
    .expect_err("prime not dividing v should be rejected");

    assert_eq!(
        error,
        EndomorphismRingLevelRecoveryError::LocalPrimeNotInFrobeniusConductor {
            prime: BigUint::from(3u8),
        }
    );
}

#[test]
fn global_report_rejects_inconsistent_local_valuation() {
    let report_for_two = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("2-local report should build");

    let error = EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set(-64),
        vec![report_for_two],
    )
    .expect_err("same prime with different v_ell(v) should be rejected");

    assert_eq!(
        error,
        EndomorphismRingLevelRecoveryError::InconsistentLocalConductorValuation {
            prime: BigUint::from(2u8),
            report_frobenius_conductor_valuation: 1,
            expected_frobenius_conductor_valuation: 2,
        }
    );
}

#[test]
fn global_report_rejects_mixed_node_reports() {
    let first = local_report_at_from_discriminant_and_path(
        BigUint::from(2u8),
        -108,
        vec![IsogenyGraphNodeId(0)],
    )
    .expect("2-local report should build");
    let second = local_report_at_from_discriminant_and_path(
        BigUint::from(3u8),
        -108,
        vec![IsogenyGraphNodeId(1)],
    )
    .expect("3-local report should build");

    let error = EndomorphismRingLevelRecoveryReport::from_local_reports(
        candidate_set(-108),
        vec![first, second],
    )
    .expect_err("reports from different nodes should be rejected");

    assert_eq!(
        error,
        EndomorphismRingLevelRecoveryError::MixedNodeReports {
            expected_node_id: IsogenyGraphNodeId(0),
            found_node_id: IsogenyGraphNodeId(1),
        }
    );
}
