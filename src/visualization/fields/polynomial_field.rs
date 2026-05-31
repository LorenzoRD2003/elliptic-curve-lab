use num_complex::Complex64;

use crate::fields::{
    complex_approx::ComplexApprox,
    polynomial_field::{PolynomialFieldElement, PolynomialModulus},
    prime_field::{Fp, FpElem},
    traits::Field,
};
use crate::polynomials::{
    DensePolynomial, IrreducibilityStatus, PolynomialError, ReducibilityReason,
    irreducibility_status,
};
use crate::visualization::Visualizable;

use super::format_complex;

/// Formats a polynomial over `GF(P)` from coefficients stored in ascending degree order.
pub fn format_prime_polynomial<const P: u64>(coefficients: &[FpElem<P>]) -> String {
    let mut terms = Vec::new();

    for (power, coefficient) in coefficients.iter().enumerate().rev() {
        if Fp::<P>::is_zero(coefficient) {
            continue;
        }

        let value = coefficient.value();
        let term = match power {
            0 => value.to_string(),
            1 if value == 1 => "x".to_string(),
            1 => format!("{value}*x"),
            _ if value == 1 => format!("x^{power}"),
            _ => format!("{value}*x^{power}"),
        };
        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

/// Formats a polynomial over the approximate complex backend.
///
/// Coefficients are interpreted in ascending degree order, exactly like the
/// dense storage used elsewhere in the crate.
pub fn format_complex_polynomial(coefficients: &[Complex64]) -> String {
    let mut terms = Vec::new();

    for (power, coefficient) in coefficients.iter().enumerate().rev() {
        if ComplexApprox::is_zero(coefficient) {
            continue;
        }

        let coeff = format_complex(coefficient);
        let term = match power {
            0 => coeff,
            1 => format!("({coeff})*x"),
            _ => format!("({coeff})*x^{power}"),
        };
        terms.push(term);
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

/// Explains how the coefficient vector maps to a polynomial over `GF(P)`.
pub fn explain_prime_polynomial_storage<const P: u64>(coefficients: &[FpElem<P>]) -> String {
    let mut lines = vec![
        format!("Polynomial over GF({P})"),
        "coefficients are stored in ascending degree order".to_string(),
        format!("polynomial: {}", format_prime_polynomial(coefficients)),
        "storage mapping:".to_string(),
    ];

    if coefficients.is_empty() {
        lines.push("- empty vector represents the zero polynomial".to_string());
        return lines.join("\n");
    }

    for (power, coefficient) in coefficients.iter().enumerate() {
        lines.push(format!(
            "- index {power}: coefficient {} multiplies x^{power}",
            coefficient.value()
        ));
    }

    lines.join("\n")
}

/// Formats a modulus polynomial used in a quotient construction over `GF(P)`.
pub fn format_prime_polynomial_modulus<const P: u64>(modulus: &PolynomialModulus<Fp<P>>) -> String {
    format!("m(x) = {}", format_prime_polynomial(modulus.coefficients()))
}

/// Returns a short textual description of a modulus polynomial over `GF(P)`.
pub fn describe_prime_polynomial_modulus<const P: u64>(
    modulus: &PolynomialModulus<Fp<P>>,
) -> String {
    format!(
        "Polynomial modulus over GF({P})\n\
         degree: {}\n\
         raw coefficients (ascending): {:?}\n\
         expression: {}",
        modulus.degree(),
        modulus
            .coefficients()
            .iter()
            .map(FpElem::value)
            .collect::<Vec<_>>(),
        format_prime_polynomial_modulus(modulus)
    )
}

/// Describes whether a modulus polynomial over `GF(P)` is suitable as the
/// defining polynomial of a quotient field.
///
/// This helper bridges the `fields` and `polynomials` educational layers:
///
/// - it starts from the field-oriented notion of a modulus polynomial
/// - it reuses the polynomial irreducibility checker currently implemented for
///   dense polynomials over prime fields
/// - it explains the practical consequence for the quotient construction
///
/// The underlying irreducibility check currently uses the library's
/// exhaustive educational baseline algorithm.
pub fn describe_prime_polynomial_modulus_as_field_modulus<const P: u64>(
    modulus: &PolynomialModulus<Fp<P>>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<Fp<P>>::new(modulus.coefficients().to_vec());
    let status = irreducibility_status(&dense_modulus)?;

    let suitability = match &status {
        IrreducibilityStatus::Irreducible => {
            "suitable for a quotient field: yes, this modulus is irreducible over the base field"
        }
        IrreducibilityStatus::Linear => {
            "suitable for a quotient field: technically yes, but the quotient is just the base field again"
        }
        IrreducibilityStatus::Reducible { .. }
        | IrreducibilityStatus::ReducibleWithoutWitness { .. } => {
            "suitable for a quotient field: no, a reducible modulus gives only a quotient algebra in general"
        }
        IrreducibilityStatus::Constant => {
            "suitable for a quotient field: no, a constant polynomial cannot define a meaningful field extension"
        }
    };

    Ok(format!(
        "Field-modulus check over GF({P})\n\
         expression: {}\n\
         base field algebraically closed: {}\n\
         irreducibility status: {}\n\
         {}",
        format_prime_polynomial_modulus(modulus),
        if Fp::<P>::IS_ALGEBRAICALLY_CLOSED {
            "yes"
        } else {
            "no"
        },
        format_irreducibility_status(&status),
        suitability
    ))
}

/// Explains the irreducibility result for a modulus polynomial over `GF(P)`.
///
/// This helper is aimed at the field-construction use case: it explains not
/// just whether the polynomial is reducible, but what that means for the
/// quotient `GF(P)[x] / (m(x))`.
pub fn explain_prime_polynomial_modulus_irreducibility<const P: u64>(
    modulus: &PolynomialModulus<Fp<P>>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<Fp<P>>::new(modulus.coefficients().to_vec());
    let status = irreducibility_status(&dense_modulus)?;

    let mut lines = vec![
        format!("Irreducibility check for a field modulus over GF({P})"),
        format!("modulus: {}", format_prime_polynomial_modulus(modulus)),
        format!(
            "base field algebraically closed: {}",
            if Fp::<P>::IS_ALGEBRAICALLY_CLOSED {
                "yes"
            } else {
                "no"
            }
        ),
    ];

    if !Fp::<P>::IS_ALGEBRAICALLY_CLOSED {
        lines.push(
            "note: the base field is not algebraically closed, so higher-degree irreducible polynomials may exist".to_string(),
        );
    }

    match status {
        IrreducibilityStatus::Constant => {
            lines.push("status: constant".to_string());
            lines.push(
                "consequence: a constant polynomial is not a valid field-extension modulus"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Linear => {
            lines.push("status: linear".to_string());
            lines.push(
                "consequence: a linear modulus is irreducible, but the quotient does not create a genuine new extension"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Irreducible => {
            lines.push("status: irreducible".to_string());
            lines.push(
                "consequence: this modulus is suitable for a quotient-field construction over the base field".to_string(),
            );
        }
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "witness divisor: {}",
                format_prime_polynomial(divisor.coefficients())
            ));
            lines.push(format!(
                "witness quotient: {}",
                format_prime_polynomial(quotient.coefficients())
            ));
            lines.push(format!(
                "factorization: {} = ({}) * ({})",
                format_prime_polynomial(modulus.coefficients()),
                format_prime_polynomial(divisor.coefficients()),
                format_prime_polynomial(quotient.coefficients())
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "current explanation: {}",
                format_reducibility_reason(reason)
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
    }

    Ok(lines.join("\n"))
}

/// Describes whether a modulus polynomial over the approximate complex
/// backend can define a non-trivial quotient field.
///
/// Because `ComplexApprox` models an algebraically closed field, every
/// polynomial of degree at least `2` is reducible. The current irreducibility
/// API reports that conclusion without fabricating a numerical factorization
/// witness.
pub fn describe_complex_polynomial_modulus_as_field_modulus(
    modulus: &PolynomialModulus<ComplexApprox>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<ComplexApprox>::new(modulus.coefficients().to_vec());
    let status = irreducibility_status(&dense_modulus)?;

    let suitability = match &status {
        IrreducibilityStatus::Irreducible => {
            "suitable for a quotient field: yes, this modulus is irreducible over the base field"
        }
        IrreducibilityStatus::Linear => {
            "suitable for a quotient field: technically yes, but the quotient is just the base field again"
        }
        IrreducibilityStatus::Reducible { .. }
        | IrreducibilityStatus::ReducibleWithoutWitness { .. } => {
            "suitable for a quotient field: no, a reducible modulus gives only a quotient algebra in general"
        }
        IrreducibilityStatus::Constant => {
            "suitable for a quotient field: no, a constant polynomial cannot define a meaningful field extension"
        }
    };

    Ok(format!(
        "Field-modulus check over C (approx)\n\
         expression: m(x) = {}\n\
         base field algebraically closed: yes\n\
         irreducibility status: {}\n\
         {}",
        format_complex_polynomial(modulus.coefficients()),
        format_complex_irreducibility_status(&status),
        suitability
    ))
}

/// Explains the irreducibility result for a modulus polynomial over the
/// approximate complex backend.
pub fn explain_complex_polynomial_modulus_irreducibility(
    modulus: &PolynomialModulus<ComplexApprox>,
) -> Result<String, PolynomialError> {
    let dense_modulus = DensePolynomial::<ComplexApprox>::new(modulus.coefficients().to_vec());
    let status = irreducibility_status(&dense_modulus)?;

    let mut lines = vec![
        "Irreducibility check for a field modulus over C (approx)".to_string(),
        format!("modulus: m(x) = {}", format_complex_polynomial(modulus.coefficients())),
        "base field algebraically closed: yes".to_string(),
        "note: in an algebraically closed field, every polynomial of degree at least two is reducible".to_string(),
    ];

    match status {
        IrreducibilityStatus::Constant => {
            lines.push("status: constant".to_string());
            lines.push(
                "consequence: a constant polynomial is not a valid field-extension modulus"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Linear => {
            lines.push("status: linear".to_string());
            lines.push(
                "consequence: a linear modulus is irreducible, but the quotient does not create a genuine new extension"
                    .to_string(),
            );
        }
        IrreducibilityStatus::Irreducible => {
            lines.push("status: irreducible".to_string());
            lines.push(
                "consequence: this modulus is suitable for a quotient-field construction over the base field".to_string(),
            );
        }
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "witness divisor: {}",
                format_complex_polynomial(divisor.coefficients())
            ));
            lines.push(format!(
                "witness quotient: {}",
                format_complex_polynomial(quotient.coefficients())
            ));
            lines.push(format!(
                "factorization: {} = ({}) * ({})",
                format_complex_polynomial(modulus.coefficients()),
                format_complex_polynomial(divisor.coefficients()),
                format_complex_polynomial(quotient.coefficients())
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            lines.push("status: reducible".to_string());
            lines.push(format!(
                "current explanation: {}",
                format_reducibility_reason(reason)
            ));
            lines.push(
                "consequence: a reducible modulus does not define a field extension in general"
                    .to_string(),
            );
        }
    }

    Ok(lines.join("\n"))
}

fn format_irreducibility_status<const P: u64>(status: &IrreducibilityStatus<Fp<P>>) -> String {
    match status {
        IrreducibilityStatus::Constant => "constant".to_string(),
        IrreducibilityStatus::Linear => "linear".to_string(),
        IrreducibilityStatus::Irreducible => "irreducible".to_string(),
        IrreducibilityStatus::Reducible { divisor, quotient } => format!(
            "reducible; witness: {} = ({}) * ({})",
            format_prime_polynomial_modulus(
                &PolynomialModulus::<Fp<P>>::new(quotient.mul(divisor).coefficients().to_vec())
                    .expect("product of non-trivial factors is a valid non-constant modulus")
            ),
            format_prime_polynomial(divisor.coefficients()),
            format_prime_polynomial(quotient.coefficients())
        ),
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            format!("reducible; reason: {}", format_reducibility_reason(*reason))
        }
    }
}

fn format_complex_irreducibility_status(status: &IrreducibilityStatus<ComplexApprox>) -> String {
    match status {
        IrreducibilityStatus::Constant => "constant".to_string(),
        IrreducibilityStatus::Linear => "linear".to_string(),
        IrreducibilityStatus::Irreducible => "irreducible".to_string(),
        IrreducibilityStatus::Reducible { divisor, quotient } => {
            let product = quotient.mul(divisor);
            format!(
                "reducible; witness: {} = ({}) * ({})",
                format_complex_polynomial(product.coefficients()),
                format_complex_polynomial(divisor.coefficients()),
                format_complex_polynomial(quotient.coefficients())
            )
        }
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
            format!("reducible; reason: {}", format_reducibility_reason(*reason))
        }
    }
}

fn format_reducibility_reason(reason: ReducibilityReason) -> &'static str {
    match reason {
        ReducibilityReason::AlgebraicallyClosedField => {
            "the base field is algebraically closed, so every degree >= 2 polynomial factors non-trivially"
        }
    }
}

/// Formats a quotient representative over `GF(P)` together with its modulus.
pub fn format_prime_polynomial_field_element<const P: u64>(
    element: &PolynomialFieldElement<Fp<P>>,
) -> String {
    format!(
        "[{}] mod ({})",
        format_prime_polynomial(element.coefficients()),
        format_prime_polynomial(element.modulus().coefficients())
    )
}

/// Returns a short educational description of a quotient element over `GF(P)`.
pub fn describe_prime_polynomial_field_element<const P: u64>(
    element: &PolynomialFieldElement<Fp<P>>,
) -> String {
    format!(
        "Quotient element over GF({P})\n\
         representative coefficients (ascending): {:?}\n\
         representative polynomial: {}\n\
         modulus polynomial: {}\n\
         note: the representative may still be unreduced in this scaffold phase",
        element
            .coefficients()
            .iter()
            .map(FpElem::value)
            .collect::<Vec<_>>(),
        format_prime_polynomial(element.coefficients()),
        format_prime_polynomial_modulus(element.modulus())
    )
}

impl<const P: u64> Visualizable for PolynomialModulus<Fp<P>> {
    fn format_compact(&self) -> String {
        format_prime_polynomial_modulus(self)
    }

    fn describe(&self) -> String {
        describe_prime_polynomial_modulus(self)
    }
}

impl<const P: u64> Visualizable for PolynomialFieldElement<Fp<P>> {
    fn format_compact(&self) -> String {
        format_prime_polynomial_field_element(self)
    }

    fn describe(&self) -> String {
        describe_prime_polynomial_field_element(self)
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::{
        describe_complex_polynomial_modulus_as_field_modulus,
        describe_prime_polynomial_field_element, describe_prime_polynomial_modulus,
        describe_prime_polynomial_modulus_as_field_modulus,
        explain_complex_polynomial_modulus_irreducibility,
        explain_prime_polynomial_modulus_irreducibility, explain_prime_polynomial_storage,
        format_complex_polynomial, format_prime_polynomial, format_prime_polynomial_field_element,
        format_prime_polynomial_modulus,
    };
    use crate::fields::{ComplexApprox, Field, Fp, PolynomialFieldElement, PolynomialModulus};
    use crate::visualization::Visualizable;

    type F17 = Fp<17>;

    fn coeffs(values: &[u64]) -> Vec<<F17 as Field>::Elem> {
        values.iter().copied().map(F17::elem_from_u64).collect()
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
        assert_eq!(format_prime_polynomial::<17>(&coeffs(&[])), "0");
        assert_eq!(format_prime_polynomial::<17>(&coeffs(&[0, 0, 0])), "0");
    }

    #[test]
    fn polynomial_formatter_handles_sparse_and_dense_terms() {
        assert_eq!(format_prime_polynomial::<17>(&coeffs(&[5])), "5");
        assert_eq!(format_prime_polynomial::<17>(&coeffs(&[0, 1])), "x");
        assert_eq!(
            format_prime_polynomial::<17>(&coeffs(&[3, 2, 0, 1])),
            "x^3 + 2*x + 3"
        );
        assert_eq!(
            format_prime_polynomial::<17>(&coeffs(&[1, 0, 4])),
            "4*x^2 + 1"
        );
    }

    #[test]
    fn polynomial_storage_explanation_mentions_order_and_indices() {
        let explanation = explain_prime_polynomial_storage::<17>(&coeffs(&[3, 0, 1]));
        assert!(explanation.contains("ascending degree order"));
        assert!(explanation.contains("index 0: coefficient 3 multiplies x^0"));
        assert!(explanation.contains("index 2: coefficient 1 multiplies x^2"));
        assert!(explanation.contains("x^2 + 3"));
    }

    #[test]
    fn modulus_visualization_is_readable() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
        assert_eq!(format_prime_polynomial_modulus(&modulus), "m(x) = x^2 + 1");

        let description = describe_prime_polynomial_modulus(&modulus);
        assert!(description.contains("Polynomial modulus over GF(17)"));
        assert!(description.contains("degree: 2"));
        assert!(description.contains("raw coefficients (ascending): [1, 0, 1]"));
        assert!(description.contains("expression: m(x) = x^2 + 1"));
    }

    #[test]
    fn quotient_element_visualization_is_readable() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
        let element = PolynomialFieldElement::<F17>::new(coeffs(&[2, 3]), modulus)
            .expect("element should exist");

        assert_eq!(
            format_prime_polynomial_field_element(&element),
            "[3*x + 2] mod (x^2 + 1)"
        );

        let description = describe_prime_polynomial_field_element(&element);
        assert!(description.contains("Quotient element over GF(17)"));
        assert!(description.contains("representative coefficients (ascending): [2, 3]"));
        assert!(description.contains("representative polynomial: 3*x + 2"));
        assert!(description.contains("modulus polynomial: m(x) = x^2 + 1"));
        assert!(description.contains("may still be unreduced"));
    }

    #[test]
    fn polynomial_visualizable_trait_reuses_core_helpers() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");
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
    }

    #[test]
    fn field_modulus_description_reports_irreducible_case() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[3, 0, 1])).expect("modulus should exist");

        let description = describe_prime_polynomial_modulus_as_field_modulus(&modulus)
            .expect("irreducibility check should work");

        assert!(description.contains("Field-modulus check over GF(17)"));
        assert!(description.contains("irreducibility status: irreducible"));
        assert!(description.contains("suitable for a quotient field: yes"));
    }

    #[test]
    fn field_modulus_description_reports_reducible_case() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");

        let description = describe_prime_polynomial_modulus_as_field_modulus(&modulus)
            .expect("irreducibility check should work");

        assert!(description.contains("irreducibility status: reducible"));
        assert!(description.contains("suitable for a quotient field: no"));
    }

    #[test]
    fn irreducibility_explanation_reports_factorization_witness() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[1, 0, 1])).expect("modulus should exist");

        let explanation = explain_prime_polynomial_modulus_irreducibility(&modulus)
            .expect("irreducibility check should work");

        assert!(explanation.contains("Irreducibility check for a field modulus over GF(17)"));
        assert!(explanation.contains("status: reducible"));
        assert!(explanation.contains("witness divisor:"));
        assert!(explanation.contains("witness quotient:"));
        assert!(
            explanation
                .contains("consequence: a reducible modulus does not define a field extension")
        );
    }

    #[test]
    fn irreducibility_explanation_reports_irreducible_case() {
        let modulus =
            PolynomialModulus::<F17>::new(coeffs(&[3, 0, 1])).expect("modulus should exist");

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
            explanation
                .contains("consequence: a reducible modulus does not define a field extension")
        );
    }
}
