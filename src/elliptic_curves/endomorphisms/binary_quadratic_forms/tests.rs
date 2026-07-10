use num_bigint::BigInt;

use crate::elliptic_curves::endomorphisms::{
    BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup,
    binary_quadratic_forms::class_group::equivalence::properly_equivalent_form,
    quadratic_orders::QuadraticDiscriminant,
};
use crate::fields::{Q, traits::Field};

fn z(value: i64) -> BigInt {
    BigInt::from(value)
}

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
fn concordant_composition_uses_the_principal_form_as_identity() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let principal = BinaryQuadraticForm::new(z(1), z(1), z(6));
    let form = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let product = class_group
        .compose_concordant_reduced_forms(&principal, &form)
        .expect("the principal form is concordant with this representative");

    assert_eq!(product, form);
}

#[test]
fn concordant_composition_reduces_the_unreduced_product() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let form = BinaryQuadraticForm::new(z(2), z(-1), z(3));

    let product = class_group
        .compose_concordant_reduced_forms(&form, &form)
        .expect("this pair satisfies gcd(a,a′,(b+b′)/2) = 1");

    assert_eq!(product, BinaryQuadraticForm::new(z(2), z(1), z(3)));
}

#[test]
fn concordant_composition_rejects_non_concordant_representatives() {
    let class_group = QuadraticClassGroup::new(QuadraticDiscriminant::new(-23))
        .expect("D = -23 should be supported");
    let left = BinaryQuadraticForm::new(z(2), z(-1), z(3));
    let right = BinaryQuadraticForm::new(z(2), z(1), z(3));

    let error = class_group
        .compose_concordant_reduced_forms(&left, &right)
        .expect_err("gcd(2,2,0) = 2, so the pair is not concordant");

    assert_eq!(error, BinaryQuadraticFormError::NotConcordantForms);
}
