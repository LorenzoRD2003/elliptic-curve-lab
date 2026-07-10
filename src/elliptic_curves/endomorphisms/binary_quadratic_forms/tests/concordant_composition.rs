use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{
        BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
        class_group::equivalence::properly_equivalent_form,
    },
    quadratic_orders::QuadraticDiscriminant,
};

use super::z;

#[test]
fn proper_equivalence_preserves_discriminant_primitivity_and_positive_definiteness() {
    let form = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let equivalent = properly_equivalent_form(&form, z(1), z(1))
        .expect("1 and 1 should complete to a proper unimodular matrix");

    assert_eq!(equivalent.discriminant(), form.discriminant());
    assert!(equivalent.is_primitive());
    assert!(equivalent.is_positive_definite());
}

#[test]
fn proper_equivalence_applies_the_expected_unimodular_substitution() {
    let form = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let equivalent = properly_equivalent_form(&form, z(1), z(1))
        .expect("1 and 1 should complete to a proper unimodular matrix");

    assert_eq!(equivalent, BinaryQuadraticForm::new(z(4), z(-3), z(2)));
}

#[test]
fn proper_equivalence_rejects_non_coprime_first_column() {
    let form = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    assert_eq!(properly_equivalent_form(&form, z(2), z(2)), None);
}

#[test]
fn public_compose_rejects_wrong_discriminant() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let left = BinaryQuadraticForm::new(z(1), z(1), z(6));
    let right = BinaryQuadraticForm::new(z(1), z(0), z(5));

    let error = class_group
        .compose(&left, &right)
        .expect_err("D = -20 does not belong to this class group");

    assert_eq!(
        error,
        BinaryQuadraticFormError::ClassGroupDiscriminantMismatch
    );
}

#[test]
fn public_compose_rejects_nonreduced_input() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let principal = BinaryQuadraticForm::new(z(1), z(1), z(6));
    let nonreduced =
        properly_equivalent_form(&BinaryQuadraticForm::new(z(2), z(-1), z(3)), z(1), z(1))
            .expect("1 and 1 should complete to a proper unimodular matrix");

    let error = class_group
        .compose(&principal, &nonreduced)
        .expect_err("the public API accepts only reduced representatives");

    assert_eq!(error, BinaryQuadraticFormError::NotReducedPositiveDefinite);
}
