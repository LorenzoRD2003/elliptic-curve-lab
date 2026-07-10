use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn order_of_reduced_form_returns_one_for_the_principal_class() {
    let class_group = class_group(-23);
    let principal = BinaryQuadraticForm::new(z(1), z(1), z(6));

    assert_eq!(
        class_group
            .order_of_reduced_form(&principal)
            .expect("principal class should have order one"),
        1
    );
}

#[test]
fn order_of_reduced_form_matches_the_gp_table_for_discriminant_minus_23() {
    let class_group = class_group(-23);
    let generator = BinaryQuadraticForm::new(z(2), z(-1), z(3));
    let inverse_generator = BinaryQuadraticForm::new(z(2), z(1), z(3));

    assert_eq!(
        class_group
            .order_of_reduced_form(&generator)
            .expect("generator should have a class order"),
        3
    );
    assert_eq!(
        class_group
            .order_of_reduced_form(&inverse_generator)
            .expect("inverse generator should have a class order"),
        3
    );
}

#[test]
fn order_of_reduced_form_validates_the_input_representative() {
    let class_group = class_group(-23);
    let wrong_discriminant = BinaryQuadraticForm::new(z(1), z(0), z(5));
    let nonreduced = BinaryQuadraticForm::new(z(4), z(-3), z(2));

    assert_eq!(
        class_group
            .order_of_reduced_form(&wrong_discriminant)
            .expect_err("order search should reject forms from another group"),
        BinaryQuadraticFormError::ClassGroupDiscriminantMismatch
    );
    assert_eq!(
        class_group
            .order_of_reduced_form(&nonreduced)
            .expect_err("order search should reject non-reduced representatives"),
        BinaryQuadraticFormError::NotReducedPositiveDefinite
    );
}

fn class_group(discriminant: i64) -> QuadraticClassGroup {
    QuadraticClassGroup::new(QuadraticDiscriminant::new(discriminant))
        .expect("test discriminant should define an imaginary quadratic order")
}
