use std::collections::BTreeMap;

use super::*;
use crate::elliptic_curves::endomorphisms::quadratic_orders::QuadraticDiscriminant;
use crate::isogenies::{
    class_group_action::{
        CraterOrientationWitness, OrientedCraterClassOrderComparisonError,
        OrientedCraterClassOrderStatus,
    },
    graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
};

fn oriented_report() -> crate::isogenies::class_group_action::OrientedLabeledCraterWalkReport {
    let graph = IsogenyGraphBuilder::new(f7_curve(), 3)
        .max_depth(2)
        .deduplicate_by_base_field_isomorphism(true)
        .build()
        .expect("small F_7 degree-three graph should build");
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
            (IsogenyGraphNodeId(1), IsogenyGraphNodeId(0)),
        ]),
    )
    .expect("certified two-node crater cycle should orient");

    labeled
        .with_user_orientation(orientation)
        .expect("orientation should attach to labeled walk")
}

#[test]
fn class_order_comparison_reports_the_toy_crater_mismatch() {
    let oriented = oriented_report();

    let comparison = oriented
        .compare_generator_order(&class_group_minus_23(), IsogenyGraphNodeId(0))
        .expect("comparison should compute");

    assert_eq!(comparison.generator_ideal().norm(), &bu(3));
    assert_eq!(comparison.generator_form(), &form(2, -1, 3));
    assert_eq!(comparison.class_order(), 3);
    assert_eq!(comparison.oriented_orbit_length(), Some(2));
    assert_eq!(
        comparison.status(),
        OrientedCraterClassOrderStatus::OrientedOrbitLengthDiffers
    );
}

#[test]
fn class_order_comparison_rejects_start_outside_the_oriented_crater() {
    let oriented = oriented_report();

    let error = oriented
        .compare_generator_order(&class_group_minus_23(), IsogenyGraphNodeId(99))
        .expect_err("start should belong to the oriented crater");

    assert_eq!(
        error,
        OrientedCraterClassOrderComparisonError::StartOutsideOrientedCrater {
            start: IsogenyGraphNodeId(99)
        }
    );
}

#[test]
fn class_order_comparison_rejects_a_mismatched_class_group() {
    let oriented = oriented_report();
    let wrong_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-20))
        .expect("D = -20 should define an imaginary quadratic class group");

    let error = oriented
        .compare_generator_order(&wrong_group, IsogenyGraphNodeId(0))
        .expect_err("class group discriminant should match the ideal order");

    assert_eq!(
        error,
        OrientedCraterClassOrderComparisonError::ClassGroupDiscriminantMismatch
    );
}
