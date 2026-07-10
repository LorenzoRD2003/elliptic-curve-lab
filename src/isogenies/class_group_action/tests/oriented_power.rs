use std::collections::BTreeMap;

use num_bigint::{BigInt, BigUint};

use super::*;
use crate::isogenies::{
    class_group_action::{
        CraterOrientationWitness, OrientedCraterPowerActionError, OrientedLabeledCraterWalkReport,
    },
    graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
};

fn oriented_report() -> OrientedLabeledCraterWalkReport {
    let graph = IsogenyGraphBuilder::new(cm_field_minus_23_curve(), 3)
        .max_depth(3)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_101 degree-three graph should build");
    let ideal = split_three_ideal();
    let crater = graph
        .volcano_crater_report(ideal.norm())
        .expect("crater report should build for the ideal norm");
    let labeled = graph
        .labeled_crater_walk_report(&class_group_minus_23(), ideal, IsogenyGraphNodeId(0))
        .expect("labeled crater walk should build");
    let orientation = CraterOrientationWitness::new(
        &crater,
        BTreeMap::from([
            (IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)),
            (IsogenyGraphNodeId(1), IsogenyGraphNodeId(2)),
            (IsogenyGraphNodeId(2), IsogenyGraphNodeId(0)),
        ]),
    )
    .expect("certified three-node crater cycle should orient");

    labeled
        .with_user_orientation(orientation)
        .expect("orientation should attach to labeled walk")
}

#[test]
fn oriented_power_zero_fixes_the_start_node() {
    let oriented = oriented_report();

    let report = oriented
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(0))
        .expect("zero exponent should apply");

    assert_eq!(report.exponent(), &BigInt::from(0));
    assert_eq!(report.path(), &[IsogenyGraphNodeId(0)]);
    assert_eq!(report.target(), IsogenyGraphNodeId(0));
}

#[test]
fn oriented_power_one_follows_the_positive_successor() {
    let oriented = oriented_report();

    let report = oriented
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(1))
        .expect("positive exponent should apply");

    assert_eq!(
        report.path(),
        &[IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]
    );
    assert_eq!(report.target(), IsogenyGraphNodeId(1));
    assert_eq!(report.generator_ideal().norm(), &BigUint::from(3u8));
    assert_eq!(report.generator_form().a(), &BigInt::from(2));
}

#[test]
fn oriented_power_negative_one_follows_the_inverse_direction() {
    let oriented = oriented_report();

    let report = oriented
        .apply_power_from(IsogenyGraphNodeId(1), BigInt::from(-1))
        .expect("negative exponent should apply");

    assert_eq!(
        report.path(),
        &[IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)]
    );
    assert_eq!(report.target(), IsogenyGraphNodeId(0));
}

#[test]
fn oriented_power_inverse_exponents_return_to_start() {
    let oriented = oriented_report();

    let forward = oriented
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(1))
        .expect("positive exponent should apply");
    let backward = oriented
        .apply_power_from(forward.target(), BigInt::from(-1))
        .expect("negative exponent should apply");

    assert_eq!(backward.target(), IsogenyGraphNodeId(0));
}

#[test]
fn oriented_power_adds_exponents_modulo_the_oriented_cycle() {
    let oriented = oriented_report();

    let first = oriented
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(1))
        .expect("positive exponent should apply");
    let second = oriented
        .apply_power_from(first.target(), BigInt::from(1))
        .expect("positive exponent should apply");
    let combined = oriented
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(2))
        .expect("combined exponent should apply");

    assert_eq!(second.target(), combined.target());
    assert_eq!(
        combined.path(),
        &[
            IsogenyGraphNodeId(0),
            IsogenyGraphNodeId(1),
            IsogenyGraphNodeId(2)
        ]
    );
}

#[test]
fn oriented_power_rejects_start_outside_the_oriented_crater() {
    let oriented = oriented_report();

    let error = oriented
        .apply_power_from(IsogenyGraphNodeId(99), BigInt::from(1))
        .expect_err("start should belong to the oriented crater");

    assert_eq!(
        error,
        OrientedCraterPowerActionError::StartOutsideOrientedCrater {
            start: IsogenyGraphNodeId(99)
        }
    );
}
