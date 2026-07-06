use num_bigint::BigUint;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::Field;
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId, endomorphisms::VolcanoSearchError,
};
use crate::numerics::PositivePrimeError;

use super::evidence::VolcanoFloorStatus;

type F11 = crate::fields::Fp11;
type F41 = crate::fields::Fp41;
type Curve11 = ShortWeierstrassCurve<F11>;
type Curve41 = ShortWeierstrassCurve<F41>;

fn f11_full_two_torsion_curve() -> Curve11 {
    Curve11::new(F11::from_i64(1), F11::from_i64(2)).expect("valid F11 curve")
}

fn f41_curve() -> Curve41 {
    Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid F41 curve")
}

fn f41_special_j_zero_curve() -> Curve41 {
    Curve41::new(F41::zero(), F41::one()).expect("valid j = 0 curve")
}

#[test]
fn floor_evidence_stays_unknown_for_unexpanded_boundary_node() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(0)
        .build()
        .expect("depth-zero graph should build");

    let evidence = graph
        .floor_evidence_at(IsogenyGraphNodeId(0), &BigUint::from(2u8))
        .expect("local evidence query should succeed");

    assert_eq!(evidence.observed_out_degree(), 0);
    assert_eq!(
        evidence.status(),
        VolcanoFloorStatus::UnknownBecausePartialGraph
    );
}

#[test]
fn floor_evidence_detects_special_j_invariants_before_degree_rules() {
    let graph = IsogenyGraphBuilder::new(f41_special_j_zero_curve(), 2)
        .max_depth(1)
        .build()
        .expect("special-j graph should build");

    let evidence = graph
        .floor_evidence_at(IsogenyGraphNodeId(0), &BigUint::from(2u8))
        .expect("local evidence query should succeed");

    assert_eq!(evidence.status(), VolcanoFloorStatus::SpecialJInvariant);
}

#[test]
fn floor_evidence_marks_expanded_low_degree_node_as_floor() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");

    let evidence = graph
        .floor_evidence_at(IsogenyGraphNodeId(0), &BigUint::from(5u8))
        .expect("local evidence query should succeed");

    assert_eq!(evidence.observed_out_degree(), 0);
    assert_eq!(evidence.status(), VolcanoFloorStatus::OnFloor);
}

#[test]
fn floor_evidence_marks_complete_ell_plus_one_degree_as_not_floor() {
    let graph = IsogenyGraphBuilder::new(f11_full_two_torsion_curve(), 2)
        .max_depth(1)
        .build()
        .expect("F11 graph with rational two-torsion should build");

    let evidence = graph
        .floor_evidence_at(IsogenyGraphNodeId(0), &BigUint::from(2u8))
        .expect("local evidence query should succeed");

    assert_eq!(evidence.observed_out_degree(), 3);
    assert_eq!(evidence.status(), VolcanoFloorStatus::NotOnFloor);
}

#[test]
fn find_floor_path_returns_the_start_when_it_already_has_floor_evidence() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");

    let report = graph
        .find_floor_path(IsogenyGraphNodeId(0), &BigUint::from(5u8))
        .expect("floor search should stop at the start");

    assert_eq!(report.prime(), &BigUint::from(5u8));
    assert_eq!(report.start(), IsogenyGraphNodeId(0));
    assert_eq!(report.floor(), IsogenyGraphNodeId(0));
    assert_eq!(report.path(), &[IsogenyGraphNodeId(0)]);
    assert_eq!(report.distance_to_floor(), 0);
}

#[test]
fn find_floor_path_refuses_to_walk_through_partial_boundary_evidence() {
    let graph = IsogenyGraphBuilder::new(f11_full_two_torsion_curve(), 2)
        .max_depth(1)
        .build()
        .expect("F11 graph with rational two-torsion should build");

    let error = graph
        .find_floor_path(IsogenyGraphNodeId(0), &BigUint::from(2u8))
        .expect_err("the deterministic path reaches an unexpanded boundary node");

    assert!(matches!(
        error,
        VolcanoSearchError::NodeNotFullyExpanded {
            node_id: IsogenyGraphNodeId(_)
        }
    ));
}

#[test]
fn floor_search_rejects_non_prime_local_parameter() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 2)
        .max_depth(1)
        .build()
        .expect("graph should build");

    let error = graph
        .floor_evidence_at(IsogenyGraphNodeId(0), &BigUint::from(4u8))
        .expect_err("composite local parameter should fail");

    assert_eq!(
        error,
        VolcanoSearchError::InvalidLocalPrime(PositivePrimeError::Composite)
    );
}
