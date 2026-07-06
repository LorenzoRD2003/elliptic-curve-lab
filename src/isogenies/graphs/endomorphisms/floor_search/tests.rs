use num_bigint::BigUint;
use proptest::prelude::*;
use proptest::strategy::ValueTree;
use proptest::test_runner::TestRunner;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::traits::Field;
use crate::isogenies::graphs::{
    IsogenyGraphBuilder, IsogenyGraphNodeId, endomorphisms::VolcanoSearchError,
};
use crate::numerics::PositivePrimeError;
use crate::proptest_support::isogenies::arb_volcanic_floor_search_case;

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
fn randomized_find_floor_path_returns_the_start_when_it_already_has_floor_evidence() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");
    let mut sampler = |_upper_bound: usize| None::<usize>;

    let report = graph
        .find_floor_path_with_sampler(IsogenyGraphNodeId(0), &BigUint::from(5u8), &mut sampler)
        .expect("floor search should stop at the start without sampling");

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
fn randomized_find_floor_path_reports_sampler_exhaustion() {
    let case = sample_non_floor_volcanic_case();
    let graph = case.graph();
    let mut sampler = |_upper_bound: usize| None::<usize>;

    let error = graph
        .find_floor_path_with_sampler(case.start(), case.prime(), &mut sampler)
        .expect_err("non-floor start should need one sampled neighbor");

    assert_eq!(
        error,
        VolcanoSearchError::SamplerExhausted {
            node_id: case.start()
        }
    );
}

#[test]
fn randomized_find_floor_path_rejects_out_of_range_sampler_indices() {
    let case = sample_non_floor_volcanic_case();
    let graph = case.graph();
    let mut sampler = |upper_bound: usize| Some(upper_bound);

    let error = graph
        .find_floor_path_with_sampler(case.start(), case.prime(), &mut sampler)
        .expect_err("sampler must honor the requested index range");

    assert!(matches!(
        error,
        VolcanoSearchError::SamplerIndexOutOfRange {
            node_id,
            sampled_index,
            upper_bound,
        } if node_id == case.start() && sampled_index == upper_bound
    ));
}

fn sample_non_floor_volcanic_case() -> crate::proptest_support::isogenies::VolcanicFloorSearchCase {
    arb_volcanic_floor_search_case()
        .prop_filter("start node should not already be on the floor", |case| {
            case.expected_distance_to_floor() > 0
        })
        .new_tree(&mut TestRunner::deterministic())
        .expect("volcanic case strategy should produce a value")
        .current()
}

fn sampler_from_indices(indices: Vec<usize>) -> impl FnMut(usize) -> Option<usize> {
    let mut indices = indices.into_iter();
    move |upper_bound: usize| indices.next().map(|index| index % upper_bound)
}

proptest! {
    #[test]
    fn deterministic_and_randomized_floor_search_reach_floor_on_generated_volcanoes(
        case in arb_volcanic_floor_search_case(),
        sampled_indices in prop::collection::vec(any::<usize>(), 8..=12),
    ) {
        let mut sampler = sampler_from_indices(sampled_indices);

        let deterministic = case
            .graph()
            .find_floor_path(case.start(), case.prime())
            .expect("deterministic walk should reach the generated volcano floor");
        let randomized = case
            .graph()
            .find_floor_path_with_sampler(case.start(), case.prime(), &mut sampler)
            .expect("randomized walk should reach the generated volcano floor");

        prop_assert_eq!(deterministic.start(), case.start());
        prop_assert_eq!(randomized.start(), case.start());
        prop_assert!(case.floor_nodes().contains(&deterministic.floor()));
        prop_assert!(case.floor_nodes().contains(&randomized.floor()));
        prop_assert!(deterministic.distance_to_floor() >= case.expected_distance_to_floor());
        prop_assert!(randomized.distance_to_floor() >= case.expected_distance_to_floor());
        prop_assert_eq!(deterministic.path().first().copied(), Some(case.start()));
        prop_assert_eq!(randomized.path().first().copied(), Some(case.start()));
        prop_assert_eq!(
            deterministic.path().len(),
            deterministic.distance_to_floor() + 1
        );
        prop_assert_eq!(randomized.path().len(), randomized.distance_to_floor() + 1);
    }

    #[test]
    fn shortest_floor_search_certifies_delta_on_generated_volcanoes(
        case in arb_volcanic_floor_search_case(),
    ) {
        let report = case
            .graph()
            .find_shortest_floor_path(case.start(), case.prime())
            .expect("shortest floor search should certify δ(v) on a complete generated volcano");

        prop_assert_eq!(report.prime(), case.prime());
        prop_assert_eq!(report.start(), case.start());
        prop_assert!(case.floor_nodes().contains(&report.floor()));
        prop_assert_eq!(report.distance_to_floor(), case.expected_distance_to_floor());
        prop_assert_eq!(report.path().first().copied(), Some(case.start()));
        prop_assert_eq!(report.path().last().copied(), Some(report.floor()));
        prop_assert_eq!(report.path().len(), report.distance_to_floor() + 1);
        prop_assert_eq!(report.level_from_total_depth(case.depth()), Some(case.start_level()));
    }
}

#[test]
fn shortest_floor_search_returns_the_start_when_it_already_has_floor_evidence() {
    let graph = IsogenyGraphBuilder::new(f41_curve(), 5)
        .max_depth(2)
        .build()
        .expect("degree-five graph with no rational kernels should build");

    let report = graph
        .find_shortest_floor_path(IsogenyGraphNodeId(0), &BigUint::from(5u8))
        .expect("shortest floor search should stop at the start");

    assert_eq!(report.prime(), &BigUint::from(5u8));
    assert_eq!(report.start(), IsogenyGraphNodeId(0));
    assert_eq!(report.floor(), IsogenyGraphNodeId(0));
    assert_eq!(report.path(), &[IsogenyGraphNodeId(0)]);
    assert_eq!(report.distance_to_floor(), 0);
    assert_eq!(report.level_from_total_depth(0), Some(0));
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
