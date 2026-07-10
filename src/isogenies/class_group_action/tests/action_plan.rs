use num_bigint::BigInt;

use crate::isogenies::class_group_action::{ClassGroupActionPlan, ClassGroupActionPlanError};

use super::{
    class_group_minus_23, class_group_minus_84, form, ramified_three_ideal_minus_84,
    split_eleven_ideal_minus_84, split_three_ideal,
};

#[test]
fn action_plan_for_identity_has_no_nonzero_factors() {
    let class_group = class_group_minus_23();
    let generator = split_three_ideal();

    let plan = ClassGroupActionPlan::from_local_ideals(&class_group, &[generator], &form(1, 1, 6))
        .expect("principal class should have a trivial plan");

    assert_eq!(plan.discriminant().value(), &BigInt::from(-23));
    assert_eq!(plan.target_form(), &form(1, 1, 6));
    assert!(plan.is_trivial());
    assert!(plan.factors().is_empty());
    assert_eq!(plan.generated_subgroup_order(), 3);
    assert_eq!(plan.ambient_class_number(), 3);
}

#[test]
fn action_plan_for_one_local_generator_has_one_factor() {
    let class_group = class_group_minus_23();
    let generator = split_three_ideal();

    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        std::slice::from_ref(&generator),
        &form(2, -1, 3),
    )
    .expect("the selected split ideal should generate its own form class");

    assert!(!plan.is_trivial());
    assert_eq!(plan.factors().len(), 1);
    assert_eq!(plan.factors()[0].ideal(), &generator);
    assert_eq!(plan.factors()[0].generator_form(), &form(2, -1, 3));
    assert_eq!(plan.factors()[0].exponent(), &BigInt::from(1));
}

#[test]
fn action_plan_factors_the_klein_product_for_discriminant_minus_84() {
    let class_group = class_group_minus_84();
    let first = split_eleven_ideal_minus_84();
    let second = ramified_three_ideal_minus_84();

    let plan = ClassGroupActionPlan::from_local_ideals(
        &class_group,
        &[first.clone(), second.clone()],
        &form(5, 4, 5),
    )
    .expect("the product class should factor through the two local generators");

    assert_eq!(plan.discriminant().value(), &BigInt::from(-84));
    assert_eq!(plan.target_form(), &form(5, 4, 5));
    assert_eq!(plan.generated_subgroup_order(), 4);
    assert_eq!(plan.ambient_class_number(), 4);
    assert_eq!(plan.factors().len(), 2);
    assert_eq!(plan.factors()[0].ideal(), &first);
    assert_eq!(plan.factors()[0].generator_form(), &form(2, 2, 11));
    assert_eq!(plan.factors()[0].exponent(), &BigInt::from(1));
    assert_eq!(plan.factors()[1].ideal(), &second);
    assert_eq!(plan.factors()[1].generator_form(), &form(3, 0, 7));
    assert_eq!(plan.factors()[1].exponent(), &BigInt::from(1));
}

#[test]
fn action_plan_rejects_target_outside_generated_subgroup() {
    let class_group = class_group_minus_84();
    let generator = split_eleven_ideal_minus_84();

    let error = ClassGroupActionPlan::from_local_ideals(&class_group, &[generator], &form(3, 0, 7))
        .expect_err("one Klein generator should not reach the other independent generator");

    assert_eq!(
        error,
        ClassGroupActionPlanError::TargetOutsideGeneratedSubgroup
    );
}

#[test]
fn action_plan_validates_local_generator_and_target_discriminants() {
    let class_group = class_group_minus_23();
    let wrong_generator = split_eleven_ideal_minus_84();

    let error =
        ClassGroupActionPlan::from_local_ideals(&class_group, &[wrong_generator], &form(1, 1, 6))
            .expect_err("local generators should live in the same discriminant");

    assert_eq!(
        error,
        ClassGroupActionPlanError::LocalGeneratorDiscriminantMismatch { index: 0 }
    );

    let generator = split_three_ideal();
    let target_error =
        ClassGroupActionPlan::from_local_ideals(&class_group, &[generator], &form(1, 0, 21))
            .expect_err("target forms should live in the same discriminant");

    assert_eq!(
        target_error,
        ClassGroupActionPlanError::TargetDiscriminantMismatch
    );
}
