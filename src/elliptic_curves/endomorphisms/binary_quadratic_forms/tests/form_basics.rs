use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::BinaryQuadraticForm;
use crate::fields::{Q, traits::Field};

use super::z;

#[test]
fn stores_the_integral_ternary_as_source_of_truth() {
    let form = BinaryQuadraticForm::new(z(2), z(3), z(5));

    assert_eq!(form.a(), &z(2));
    assert_eq!(form.b(), &z(3));
    assert_eq!(form.c(), &z(5));
    assert_eq!(form.coefficients(), (&z(2), &z(3), &z(5)));
}

#[test]
fn polynomial_view_is_a_binary_quadratic_over_q() {
    let form = BinaryQuadraticForm::new(z(2), z(3), z(5));
    let polynomial = form.polynomial();

    assert_eq!(polynomial.arity(), 2);
    assert_eq!(polynomial.degree(), Some(2));
    assert_eq!(polynomial.len(), 3);

    let value = polynomial
        .evaluate(&[Q::from_bigint(&z(7)), Q::from_bigint(&z(11))])
        .expect("arity should match");

    assert_eq!(
        value,
        Q::from_bigint(&form.evaluate_integral(&z(7), &z(11)))
    );
}

#[test]
fn polynomial_view_drops_zero_terms_but_keeps_arity() {
    let form = BinaryQuadraticForm::new(z(1), z(0), z(0));
    let polynomial = form.polynomial();

    assert_eq!(polynomial.arity(), 2);
    assert_eq!(polynomial.len(), 1);
    assert_eq!(polynomial.degree(), Some(2));
}

#[test]
fn discriminant_primitivity_and_positive_definiteness_are_derived_from_coefficients() {
    let principal = BinaryQuadraticForm::new(z(1), z(1), z(6));
    let imprimitive = BinaryQuadraticForm::new(z(2), z(2), z(12));
    let indefinite = BinaryQuadraticForm::new(z(1), z(0), z(-1));

    assert_eq!(principal.discriminant(), z(-23));
    assert!(principal.is_primitive());
    assert!(principal.is_positive_definite());

    assert_eq!(imprimitive.discriminant(), z(-92));
    assert!(!imprimitive.is_primitive());
    assert!(imprimitive.is_positive_definite());

    assert_eq!(indefinite.discriminant(), z(4));
    assert!(indefinite.is_primitive());
    assert!(!indefinite.is_positive_definite());
}

#[test]
fn zero_form_is_not_primitive_or_positive_definite() {
    let zero = BinaryQuadraticForm::new(z(0), z(0), z(0));

    assert_eq!(zero.discriminant(), z(0));
    assert!(!zero.is_primitive());
    assert!(!zero.is_positive_definite());
}

#[test]
fn conjugation_negates_the_middle_coefficient() {
    let form = BinaryQuadraticForm::new(z(3), z(5), z(7));
    let conjugate = form.conjugate();

    assert_eq!(conjugate.coefficients(), (&z(3), &z(-5), &z(7)));
    assert_eq!(conjugate.discriminant(), form.discriminant());
    assert_eq!(conjugate.conjugate(), form);
}

#[test]
fn evaluates_integral_coordinates_exactly() {
    let form = BinaryQuadraticForm::new(z(2), z(-3), z(5));

    assert_eq!(form.evaluate_integral(&z(4), &z(-2)), z(76));
}
