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

fn class_group(discriminant: i64) -> QuadraticClassGroup {
    QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))
        .expect("test discriminant should define an imaginary quadratic order")
}
