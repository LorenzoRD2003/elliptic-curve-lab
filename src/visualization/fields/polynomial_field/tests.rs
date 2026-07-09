use num_complex::Complex64;

use super::*;
use crate::fields::{
    complex_approx::ComplexApprox,
    polynomial_field::{PolynomialFieldElement, PolynomialModulus},
    traits::Field,
};
use crate::visualization::{Visualizable, VisualizableField};

type F17 = crate::fields::Fp17;

fn coeffs(values: &[u64]) -> Vec<<F17 as Field>::Elem> {
    values.iter().copied().map(F17::from_i64).collect()
}

fn complex_coeffs(values: &[(f64, f64)]) -> Vec<Complex64> {
    values
        .iter()
        .copied()
        .map(|(re, im)| Complex64::new(re, im))
        .collect()
}

#[test]
fn polynomial_formatter_handles_zero_polynomial() {
    assert_eq!(format_prime_polynomial::<F17>(&coeffs(&[])), "0");
    assert_eq!(format_prime_polynomial::<F17>(&coeffs(&[0, 0, 0])), "0");
}

#[test]
fn polynomial_formatter_handles_sparse_and_dense_terms() {
    assert_eq!(format_prime_polynomial::<F17>(&coeffs(&[5])), "5");
    assert_eq!(format_prime_polynomial::<F17>(&coeffs(&[0, 1])), "x");
    assert_eq!(
        format_prime_polynomial::<F17>(&coeffs(&[3, 2, 0, 1])),
        "x^3 + 2*x + 3"
    );
    assert_eq!(
        format_prime_polynomial::<F17>(&coeffs(&[1, 0, 4])),
        "4*x^2 + 1"
    );
}

#[test]
fn polynomial_storage_explanation_mentions_order_and_indices() {
    let explanation = explain_prime_polynomial_storage::<F17>(&coeffs(&[3, 0, 1]));
    assert!(explanation.contains("ascending degree order"));
    assert!(explanation.contains("index 0: coefficient 3 multiplies x^0"));
    assert!(explanation.contains("index 2: coefficient 1 multiplies x^2"));
    assert!(explanation.contains("x^2 + 3"));
}

#[test]
fn modulus_visualization_is_readable() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
    assert_eq!(format_prime_polynomial_modulus(&modulus), "m(x) = x^2 + 1");

    let description = describe_prime_polynomial_modulus(&modulus);
    assert!(description.contains("Polynomial modulus over GF(17)"));
    assert!(description.contains("degree: 2"));
    assert!(description.contains("raw coefficients (ascending): [1, 0, 1]"));
    assert!(description.contains("expression: m(x) = x^2 + 1"));
}

#[test]
fn quotient_element_visualization_is_readable() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
    let element =
        PolynomialFieldElement::<F17>::new(coeffs(&[2, 3]), modulus).expect("element should exist");

    assert_eq!(
        format_prime_polynomial_field_element(&element),
        "[3*x + 2] mod (x^2 + 1)"
    );

    let description = describe_prime_polynomial_field_element(&element);
    assert!(description.contains("Quotient element over GF(17)"));
    assert!(description.contains("representative coefficients (ascending): [2, 3]"));
    assert!(description.contains("representative polynomial: 3*x + 2"));
    assert!(description.contains("reduced representative: 3*x + 2"));
    assert!(description.contains("already reduced: yes"));
    assert!(description.contains("modulus polynomial: m(x) = x^2 + 1"));
    assert!(description.contains("arithmetic is interpreted modulo"));
}

#[test]
fn polynomial_visualizable_trait_reuses_core_helpers() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
    let element = PolynomialFieldElement::<F17>::new(coeffs(&[2, 3]), modulus.clone())
        .expect("element should exist");

    assert_eq!(modulus.format_compact(), "m(x) = x^2 + 1");
    assert!(
        modulus
            .describe()
            .contains("Polynomial modulus over GF(17)")
    );

    assert_eq!(element.format_compact(), "[3*x + 2] mod (x^2 + 1)");
    assert!(element.describe().contains("Quotient element over GF(17)"));
    assert!(VisualizableField::inverse(&element).is_some());
}

#[test]
fn quotient_reduction_explanation_reports_raw_and_reduced_forms() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
    let element = PolynomialFieldElement::<F17>::new(coeffs(&[1, 2, 3]), modulus)
        .expect("element should exist");

    let explanation = explain_prime_polynomial_field_reduction(&element)
        .expect("reduction explanation should succeed");

    assert!(explanation.contains("Reduction in GF(17)[x] / (m(x))"));
    assert!(explanation.contains("raw representative: 3*x^2 + 2*x + 1"));
    assert!(explanation.contains("reduced representative: 2*x + 15"));
}

#[test]
fn quotient_addition_and_multiplication_explanations_show_reduced_results() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
    let left = PolynomialFieldElement::<F17>::new(coeffs(&[1, 1]), modulus.clone())
        .expect("left should exist");
    let right =
        PolynomialFieldElement::<F17>::new(coeffs(&[3, 16]), modulus).expect("right should exist");

    let add = explain_prime_polynomial_field_add(&left, &right)
        .expect("addition explanation should succeed");
    let mul = explain_prime_polynomial_field_mul(&left, &right)
        .expect("multiplication explanation should succeed");

    assert!(add.contains("Addition in GF(17)[x] / (m(x))"));
    assert!(add.contains("reduced result: [4] mod (x^2 + 1)"));

    assert!(mul.contains("Multiplication in GF(17)[x] / (m(x))"));
    assert!(mul.contains("reduced result: [2*x + 4] mod (x^2 + 1)"));
}

#[test]
fn quotient_inverse_explanation_and_visualizable_trait_work() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[3, 0, 1])).expect("modulus should exist");
    let element =
        PolynomialFieldElement::<F17>::new(coeffs(&[1, 1]), modulus).expect("element should exist");

    let explanation = explain_prime_polynomial_field_inverse(&element)
        .expect("inverse explanation should succeed");
    assert!(explanation.contains("Inverse in GF(17)[x] / (m(x))"));
    assert!(explanation.contains("verification:"));

    let add = PolynomialFieldElement::<F17>::explain_add(&element, &element)
        .expect("visualizable addition should exist");
    let mul = PolynomialFieldElement::<F17>::explain_mul(&element, &element)
        .expect("visualizable multiplication should exist");
    let div = PolynomialFieldElement::<F17>::explain_div(&element, &element)
        .expect("visualizable division should exist");

    assert!(add.contains("Addition in GF(17)[x] / (m(x))"));
    assert!(mul.contains("Multiplication in GF(17)[x] / (m(x))"));
    assert!(div.contains("Division in GF(17)[x] / (m(x))"));
}

#[test]
fn field_modulus_description_reports_irreducible_case() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[3, 0, 1])).expect("modulus should exist");

    let description = describe_prime_polynomial_modulus_as_field_modulus(&modulus)
        .expect("irreducibility check should work");

    assert!(description.contains("Field-modulus check over GF(17)"));
    assert!(description.contains("irreducibility status: irreducible"));
    assert!(description.contains("suitable for a quotient field: yes"));
}

#[test]
fn field_modulus_description_reports_reducible_case() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");

    let description = describe_prime_polynomial_modulus_as_field_modulus(&modulus)
        .expect("irreducibility check should work");

    assert!(description.contains("irreducibility status: reducible"));
    assert!(description.contains("suitable for a quotient field: no"));
}

#[test]
fn irreducibility_explanation_reports_factorization_witness() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");

    let explanation = explain_prime_polynomial_modulus_irreducibility(&modulus)
        .expect("irreducibility check should work");

    assert!(explanation.contains("Irreducibility check for a field modulus over GF(17)"));
    assert!(explanation.contains("status: reducible"));
    assert!(explanation.contains("witness divisor:"));
    assert!(explanation.contains("witness quotient:"));
    assert!(
        explanation.contains("consequence: a reducible modulus does not define a field extension")
    );
}

#[test]
fn irreducibility_explanation_reports_irreducible_case() {
    let modulus = PolynomialModulus::<F17>::new(coeffs(&[3, 0, 1])).expect("modulus should exist");

    let explanation = explain_prime_polynomial_modulus_irreducibility(&modulus)
        .expect("irreducibility check should work");

    assert!(explanation.contains("status: irreducible"));
    assert!(explanation.contains("suitable for a quotient-field construction"));
}

#[test]
fn complex_polynomial_formatter_is_readable() {
    let formatted = format_complex_polynomial(&complex_coeffs(&[(1.0, 0.0), (0.0, 1.0)]));
    assert_eq!(formatted, "(0.000000 + 1.000000i)*x + 1.000000 + 0.000000i");
}

#[test]
fn complex_field_modulus_description_reports_reducibility_by_field_property() {
    let modulus = PolynomialModulus::<ComplexApprox>::new(complex_coeffs(&[
        (1.0, 0.0),
        (0.0, 0.0),
        (1.0, 0.0),
    ]))
    .expect("modulus should exist");

    let description = describe_complex_polynomial_modulus_as_field_modulus(&modulus)
        .expect("complex irreducibility check should work");

    assert!(description.contains("Field-modulus check over C (approx)"));
    assert!(description.contains("base field algebraically closed: yes"));
    assert!(description.contains("irreducibility status: reducible"));
    assert!(description.contains("algebraically closed"));
    assert!(description.contains("suitable for a quotient field: no"));
}

#[test]
fn complex_irreducibility_explanation_reports_theoretical_reducibility() {
    let modulus = PolynomialModulus::<ComplexApprox>::new(complex_coeffs(&[
        (1.0, 0.0),
        (0.0, 0.0),
        (1.0, 0.0),
    ]))
    .expect("modulus should exist");

    let explanation = explain_complex_polynomial_modulus_irreducibility(&modulus)
        .expect("complex irreducibility check should work");

    assert!(explanation.contains("Irreducibility check for a field modulus over C (approx)"));
    assert!(explanation.contains("status: reducible"));
    assert!(explanation.contains("every polynomial of degree at least two is reducible"));
    assert!(explanation.contains("current explanation:"));
    assert!(
        explanation.contains("consequence: a reducible modulus does not define a field extension")
    );
}
