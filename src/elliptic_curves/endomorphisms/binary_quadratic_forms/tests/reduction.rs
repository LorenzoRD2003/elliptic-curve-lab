use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::{
    BinaryQuadraticForm, BinaryQuadraticFormError,
};

use super::z;

#[test]
fn reduced_positive_definite_predicate_uses_gauss_boundary_convention() {
    let reduced = BinaryQuadraticForm::new(z(1), z(1), z(6));
    let negative_boundary = BinaryQuadraticForm::new(z(2), z(-2), z(3));
    let equal_sides_negative_middle = BinaryQuadraticForm::new(z(3), z(-1), z(3));
    let indefinite = BinaryQuadraticForm::new(z(1), z(0), z(-1));

    assert!(reduced.is_reduced_positive_definite());
    assert!(!negative_boundary.is_reduced_positive_definite());
    assert!(!equal_sides_negative_middle.is_reduced_positive_definite());
    assert!(!indefinite.is_reduced_positive_definite());
}

#[test]
fn reduce_positive_definite_reduces_unbalanced_forms() {
    let form = BinaryQuadraticForm::new(z(5), z(7), z(3));
    let reduced = form
        .reduce_positive_definite()
        .expect("positive-definite form should reduce");

    assert_eq!(reduced.coefficients(), (&z(1), &z(1), &z(3)));
    assert_eq!(reduced.discriminant(), form.discriminant());
    assert!(reduced.is_reduced_positive_definite());
}

#[test]
fn reduce_positive_definite_handles_boundary_normalization() {
    let negative_boundary = BinaryQuadraticForm::new(z(2), z(-2), z(3));
    let equal_sides_negative_middle = BinaryQuadraticForm::new(z(3), z(-1), z(3));

    assert_eq!(
        negative_boundary
            .reduce_positive_definite()
            .expect("positive-definite form should reduce")
            .coefficients(),
        (&z(2), &z(2), &z(3))
    );
    assert_eq!(
        equal_sides_negative_middle
            .reduce_positive_definite()
            .expect("positive-definite form should reduce")
            .coefficients(),
        (&z(3), &z(1), &z(3))
    );
}

#[test]
fn reduce_positive_definite_rejects_non_positive_definite_forms() {
    let error = BinaryQuadraticForm::new(z(1), z(0), z(-1))
        .reduce_positive_definite()
        .expect_err("indefinite form should not reduce as positive definite");

    assert_eq!(error, BinaryQuadraticFormError::NotPositiveDefinite);
}
