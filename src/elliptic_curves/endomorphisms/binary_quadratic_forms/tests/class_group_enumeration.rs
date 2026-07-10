use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn class_group_rejects_non_imaginary_or_non_order_discriminants() {
    let non_negative = QuadraticClassGroup::new(QuadraticDiscriminant::new(0))
        .expect_err("D = 0 should not define an imaginary class group");
    let bad_congruence = QuadraticClassGroup::new(QuadraticDiscriminant::new(-2))
        .expect_err("D = -2 is not congruent to 0 or 1 modulo 4");

    assert_eq!(
        non_negative,
        BinaryQuadraticFormError::NotNegativeDiscriminant
    );
    assert_eq!(
        bad_congruence,
        BinaryQuadraticFormError::NotQuadraticOrderDiscriminant
    );
}

#[test]
fn class_group_enumerates_reduced_forms_for_class_number_one_examples() {
    let minus_three = QuadraticClassGroup::new(QuadraticDiscriminant::new(-3))
        .expect("D = -3 should be supported");
    let minus_four = QuadraticClassGroup::new(QuadraticDiscriminant::new(-4))
        .expect("D = -4 should be supported");

    assert_eq!(
        minus_three.enumerate_reduced_forms(),
        vec![BinaryQuadraticForm::new(z(1), z(1), z(1))]
    );
    assert_eq!(
        minus_four.enumerate_reduced_forms(),
        vec![BinaryQuadraticForm::new(z(1), z(0), z(1))]
    );
}

#[test]
fn class_group_enumerates_primitive_reduced_forms_for_minus_twenty_three() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let forms = class_group.enumerate_reduced_forms();

    assert_eq!(
        forms,
        vec![
            BinaryQuadraticForm::new(z(1), z(1), z(6)),
            BinaryQuadraticForm::new(z(2), z(-1), z(3)),
            BinaryQuadraticForm::new(z(2), z(1), z(3)),
        ]
    );
    assert!(forms.iter().all(BinaryQuadraticForm::is_primitive));
    assert!(
        forms
            .iter()
            .all(BinaryQuadraticForm::is_reduced_positive_definite)
    );
    assert!(forms.iter().all(|form| form.discriminant() == z(-23)));
}

#[test]
fn class_group_enumerates_nonfundamental_order_discriminants() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-20))
        .expect("D = -20 should be supported");

    assert_eq!(
        class_group.enumerate_reduced_forms(),
        vec![
            BinaryQuadraticForm::new(z(1), z(0), z(5)),
            BinaryQuadraticForm::new(z(2), z(2), z(3)),
        ]
    );
}
