use crate::fields::{
    polynomial_field::{PolynomialFieldElement, PolynomialModulus},
    prime_field::{Fp, FpElem},
    traits::Field,
};
use crate::visualization::Visualizable;

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
    use super::{
        describe_prime_polynomial_field_element, describe_prime_polynomial_modulus,
        explain_prime_polynomial_storage, format_prime_polynomial,
        format_prime_polynomial_field_element, format_prime_polynomial_modulus,
    };
    use crate::fields::{Field, Fp, PolynomialFieldElement, PolynomialModulus};
    use crate::visualization::Visualizable;

    type F17 = Fp<17>;

    fn coeffs(values: &[u64]) -> Vec<<F17 as Field>::Elem> {
        values.iter().copied().map(F17::elem_from_u64).collect()
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
}
