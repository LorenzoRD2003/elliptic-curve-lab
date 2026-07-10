use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn generated_subgroup_records_powers_for_the_discriminant_minus_23_generator() {
    let class_group = class_group(-23);
    let generator = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let subgroup = class_group
        .generated_subgroup(&generator)
        .expect("D = -23 generator should produce a cyclic subgroup");

    assert_eq!(subgroup.discriminant().value(), &z(-23));
    assert_eq!(subgroup.generator(), &generator);
    assert_eq!(
        subgroup.elements(),
        &[
            BinaryQuadraticForm::new(z(1), z(1), z(6)),
            BinaryQuadraticForm::new(z(2), z(-1), z(3)),
            BinaryQuadraticForm::new(z(2), z(1), z(3)),
        ]
    );
    assert_eq!(subgroup.order(), 3);
    assert_eq!(subgroup.class_number(), 3);
    assert!(subgroup.is_whole_class_group());
    assert!(subgroup.contains_reduced_form(&BinaryQuadraticForm::new(z(2), z(1), z(3))));
}

#[test]
fn generated_subgroup_can_be_proper_inside_a_non_cyclic_class_group() {
    let class_group = class_group(-84);
    let generator = BinaryQuadraticForm::new(z(2), z(2), z(11));

    let subgroup = class_group
        .generated_subgroup(&generator)
        .expect("an order-two class should generate a proper subgroup");

    assert_eq!(
        subgroup.elements(),
        &[
            BinaryQuadraticForm::new(z(1), z(0), z(21)),
            BinaryQuadraticForm::new(z(2), z(2), z(11)),
        ]
    );
    assert_eq!(subgroup.order(), 2);
    assert_eq!(subgroup.class_number(), 4);
    assert!(!subgroup.is_whole_class_group());
    assert!(!subgroup.contains_reduced_form(&BinaryQuadraticForm::new(z(3), z(0), z(7))));
}

#[test]
fn generated_subgroup_validates_the_generator() {
    let class_group = class_group(-23);
    let wrong_discriminant = BinaryQuadraticForm::new(z(1), z(0), z(5));
    let nonreduced = BinaryQuadraticForm::new(z(4), z(-3), z(2));

    assert_eq!(
        class_group
            .generated_subgroup(&wrong_discriminant)
            .expect_err("subgroup generation should reject forms from another group"),
        BinaryQuadraticFormError::ClassGroupDiscriminantMismatch
    );
    assert_eq!(
        class_group
            .generated_subgroup(&nonreduced)
            .expect_err("subgroup generation should reject non-reduced representatives"),
        BinaryQuadraticFormError::NotReducedPositiveDefinite
    );
}

#[test]
fn generated_subgroup_by_set_matches_the_cyclic_report_for_one_generator() {
    let class_group = class_group(-23);
    let generator = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let cyclic = class_group
        .generated_subgroup(&generator)
        .expect("single generator should produce a cyclic report");
    let generated = class_group
        .generated_subgroup_by_set(std::slice::from_ref(&generator))
        .expect("single generator should also produce a set-generated report");

    assert_eq!(generated.discriminant(), cyclic.discriminant());
    assert_eq!(generated.generators(), std::slice::from_ref(&generator));
    assert_eq!(generated.elements(), cyclic.elements());
    assert_eq!(generated.order(), cyclic.order());
    assert_eq!(generated.class_number(), cyclic.class_number());
    assert!(generated.is_whole_class_group());
}

#[test]
fn generated_subgroup_by_set_recovers_the_klein_group_for_discriminant_minus_84() {
    let class_group = class_group(-84);
    let first = BinaryQuadraticForm::new(z(2), z(2), z(11));
    let second = BinaryQuadraticForm::new(z(3), z(0), z(7));

    let generated = class_group
        .generated_subgroup_by_set(&[first.clone(), second.clone()])
        .expect("two independent order-two classes should generate the Klein group");

    assert_eq!(generated.discriminant().value(), &z(-84));
    assert_eq!(generated.generators(), &[first, second]);
    assert_eq!(
        generated.elements(),
        &[
            BinaryQuadraticForm::new(z(1), z(0), z(21)),
            BinaryQuadraticForm::new(z(2), z(2), z(11)),
            BinaryQuadraticForm::new(z(3), z(0), z(7)),
            BinaryQuadraticForm::new(z(5), z(4), z(5)),
        ]
    );
    assert_eq!(generated.order(), 4);
    assert_eq!(generated.class_number(), 4);
    assert!(generated.is_whole_class_group());
    assert!(generated.contains_reduced_form(&BinaryQuadraticForm::new(z(5), z(4), z(5))));
}

#[test]
fn generated_subgroup_by_set_can_still_be_proper_for_one_klein_generator() {
    let class_group = class_group(-84);
    let generator = BinaryQuadraticForm::new(z(2), z(2), z(11));

    let generated = class_group
        .generated_subgroup_by_set(&[generator])
        .expect("one order-two generator should produce a proper subgroup");

    assert_eq!(generated.order(), 2);
    assert_eq!(generated.class_number(), 4);
    assert!(!generated.is_whole_class_group());
    assert!(!generated.contains_reduced_form(&BinaryQuadraticForm::new(z(3), z(0), z(7))));
}

#[test]
fn generated_subgroup_by_set_validates_the_generators() {
    let class_group = class_group(-23);
    let generator = BinaryQuadraticForm::new(z(2), z(-1), z(3));
    let wrong_discriminant = BinaryQuadraticForm::new(z(1), z(0), z(5));
    let nonreduced = BinaryQuadraticForm::new(z(4), z(-3), z(2));

    assert_eq!(
        class_group
            .generated_subgroup_by_set(&[])
            .expect_err("empty generator sets should be rejected"),
        BinaryQuadraticFormError::EmptyGeneratorSet
    );
    assert_eq!(
        class_group
            .generated_subgroup_by_set(&[generator.clone(), wrong_discriminant])
            .expect_err("set generation should reject forms from another group"),
        BinaryQuadraticFormError::ClassGroupDiscriminantMismatch
    );
    assert_eq!(
        class_group
            .generated_subgroup_by_set(&[generator, nonreduced])
            .expect_err("set generation should reject non-reduced representatives"),
        BinaryQuadraticFormError::NotReducedPositiveDefinite
    );
}

fn class_group(discriminant: i64) -> QuadraticClassGroup {
    QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))
        .expect("test discriminant should define an imaginary quadratic order")
}
