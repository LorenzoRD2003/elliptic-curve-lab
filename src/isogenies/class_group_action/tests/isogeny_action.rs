use num_bigint::BigInt;
use proptest::prelude::*;

use super::*;
use crate::isogenies::{
    class_group_action::{
        ClassGroupActionPlan, ClassGroupIsogenyActionError, CraterDirectionCertification,
        CraterOrientationWitness, OrientedCraterPowerActionError, OrientedLabeledCraterWalkReport,
    },
    graphs::{IsogenyGraphBuilder, IsogenyGraphNodeId},
};

fn oriented_split_three_report() -> OrientedLabeledCraterWalkReport {
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
    let orientation = CraterOrientationWitness::new(&crater, orientation_012())
        .expect("certified three-node crater cycle should orient");

    labeled
        .with_user_orientation(orientation)
        .expect("orientation should attach to labeled walk")
}

#[test]
fn trivial_action_plan_execution_stays_at_the_start() {
    let class_group = class_group_minus_23();
    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[split_three_ideal()],
        &form(1, 1, 6),
    )
    .expect("principal class should have a trivial plan");

    let report = plan
        .execute_from(IsogenyGraphNodeId(2), &[])
        .expect("trivial plan should need no witnesses");

    assert!(report.is_trivial());
    assert_eq!(report.start(), IsogenyGraphNodeId(2));
    assert_eq!(report.target(), IsogenyGraphNodeId(2));
    assert!(report.segments().is_empty());
}

#[test]
fn action_plan_execution_applies_the_matched_oriented_local_factor() {
    let class_group = class_group_minus_23();
    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[split_three_ideal()],
        &form(2, -1, 3),
    )
    .expect("split ideal should produce a one-factor plan");
    let witness = oriented_split_three_report();

    let report = plan
        .execute_from(IsogenyGraphNodeId(0), &[witness])
        .expect("matching oriented local witness should execute the plan");

    assert_eq!(report.start(), IsogenyGraphNodeId(0));
    assert_eq!(report.target(), IsogenyGraphNodeId(1));
    assert_eq!(report.segments().len(), 1);
    let segment = &report.segments()[0];
    assert_eq!(segment.factor_index(), 0);
    assert_eq!(segment.ideal(), &split_three_ideal());
    assert_eq!(segment.generator_form(), &form(2, -1, 3));
    assert_eq!(segment.exponent(), &BigInt::from(1));
    assert_eq!(
        segment.path(),
        &[IsogenyGraphNodeId(0), IsogenyGraphNodeId(1)]
    );
    assert_eq!(segment.start(), IsogenyGraphNodeId(0));
    assert_eq!(segment.target(), IsogenyGraphNodeId(1));
    assert_eq!(
        segment.direction_certification(),
        CraterDirectionCertification::UserSuppliedArithmeticOrientation
    );
}

#[test]
fn action_plan_execution_rejects_missing_local_witness() {
    let class_group = class_group_minus_23();
    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[split_three_ideal()],
        &form(2, -1, 3),
    )
    .expect("split ideal should produce a one-factor plan");

    let error = plan
        .execute_from(IsogenyGraphNodeId(0), &[])
        .expect_err("nontrivial plan needs an oriented local witness");

    assert_eq!(
        error,
        ClassGroupIsogenyActionError::MissingLocalWitness {
            factor_index: 0,
            ideal_norm: bu(3),
            generator_form: form(2, -1, 3)
        }
    );
}

#[test]
fn action_plan_execution_rejects_start_outside_matched_oriented_crater() {
    let class_group = class_group_minus_23();
    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[split_three_ideal()],
        &form(2, -1, 3),
    )
    .expect("split ideal should produce a one-factor plan");
    let witness = oriented_split_three_report();

    let error = plan
        .execute_from(IsogenyGraphNodeId(99), &[witness])
        .expect_err("execution should stay inside each oriented crater");

    assert_eq!(
        error,
        ClassGroupIsogenyActionError::LocalPower {
            factor_index: 0,
            source: OrientedCraterPowerActionError::StartOutsideOrientedCrater {
                start: IsogenyGraphNodeId(99)
            }
        }
    );
}

#[test]
fn executed_local_factor_matches_direct_oriented_power() {
    let class_group = class_group_minus_23();
    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[split_three_ideal()],
        &form(2, 1, 3),
    )
    .expect("square of the split ideal should produce exponent two");
    let witness = oriented_split_three_report();
    let direct = witness
        .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(2))
        .expect("direct local power should apply");

    let report = plan
        .execute_from(IsogenyGraphNodeId(0), &[witness])
        .expect("matching oriented local witness should execute the plan");

    assert_eq!(report.target(), direct.target());
    assert_eq!(report.segments()[0].path(), direct.path());
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(12))]

    #[test]
    fn executed_local_factor_agrees_with_direct_power_for_small_exponents(exponent in 0usize..9) {
        let class_group = class_group_minus_23();
        let target = match exponent % 3 {
            0 => form(1, 1, 6),
            1 => form(2, -1, 3),
            _ => form(2, 1, 3),
        };
        let plan = ClassGroupActionPlan::from_local_ideals(
            &class_group,
            &[split_three_ideal()],
            &target,
        )
        .expect("cyclic D = -23 target should factor through the split ideal");
        let witness = oriented_split_three_report();
        let direct = witness
            .apply_power_from(IsogenyGraphNodeId(0), BigInt::from(exponent))
            .expect("direct local power should apply");

        let report = plan
            .execute_from(IsogenyGraphNodeId(0), &[witness])
            .expect("matching oriented local witness should execute the plan");

        prop_assert_eq!(report.target(), direct.target());
        if let Some(segment) = report.segments().first() {
            prop_assert_eq!(segment.path(), direct.path());
        } else {
            prop_assert_eq!(direct.path(), &[IsogenyGraphNodeId(0)]);
        }
    }
}
