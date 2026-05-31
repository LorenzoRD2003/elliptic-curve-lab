use crate::fields::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};
use crate::visualization::VisualizableField;

use super::format_dense_polynomial;

/// Explains univariate Euclidean division for dense polynomials over a field.
///
/// The explanation follows the classical long-division process step by step:
///
/// - inspect the current remainder
/// - match its leading term against the divisor's leading term
/// - form the next quotient term
/// - subtract the appropriate shifted multiple of the divisor
///
/// At the end, the explanation reports the final quotient and remainder and
/// checks the recomposition identity
///
/// `dividend = divisor * quotient + remainder`.
pub fn explain_dense_division<F>(
    dividend: &DensePolynomial<F>,
    divisor: &DensePolynomial<F>,
) -> Result<String, PolynomialError>
where
    F: Field,
    F::Elem: VisualizableField,
{
    if divisor.is_zero() {
        return Err(PolynomialError::DivisionByZeroPolynomial);
    }

    let mut lines = vec![
        "Dense polynomial division".to_string(),
        format!("dividend: {}", format_dense_polynomial(dividend)),
        format!("divisor: {}", format_dense_polynomial(divisor)),
    ];

    let (quotient, remainder) = dividend.div_rem(divisor)?;

    if dividend.is_zero() {
        lines.push("the dividend is zero, so quotient and remainder are both zero".to_string());
    } else if dividend.degree() < divisor.degree() {
        lines.push(
            "the divisor degree is larger than the dividend degree, so the quotient is zero and the remainder is the dividend"
                .to_string(),
        );
    } else {
        let divisor_degree = divisor.degree().expect("non-zero divisor has degree");
        let divisor_leading = divisor
            .leading_coefficient()
            .expect("non-zero divisor has leading coefficient")
            .clone();

        let mut running_quotient = DensePolynomial::<F>::new(Vec::new());
        let mut running_remainder = dividend.clone();
        let mut step = 0usize;

        while let Some(remainder_degree) = running_remainder.degree() {
            if remainder_degree < divisor_degree {
                break;
            }

            let remainder_leading = running_remainder
                .leading_coefficient()
                .expect("non-zero remainder has leading coefficient")
                .clone();
            let degree_gap = remainder_degree - divisor_degree;
            let scale = F::div(&remainder_leading, &divisor_leading)
                .map_err(|_| PolynomialError::NonInvertibleLeadingCoefficient)?;

            let quotient_term = build_monomial::<F>(scale.clone(), degree_gap);
            let subtraction_term = shift_dense(&divisor.scale(&scale), degree_gap);
            let next_remainder = running_remainder.sub(&subtraction_term);
            running_quotient = running_quotient.add(&quotient_term);

            lines.push(format!("step {step}:"));
            lines.push(format!(
                "- current remainder: {}",
                format_dense_polynomial(&running_remainder)
            ));
            lines.push(format!(
                "- leading coefficient ratio: {} / {} = {}",
                remainder_leading.format_elem(),
                divisor_leading.format_elem(),
                scale.format_elem()
            ));
            lines.push(format!(
                "- degree gap: {} - {} = {}",
                remainder_degree, divisor_degree, degree_gap
            ));
            lines.push(format!(
                "- next quotient term: {}",
                format_dense_polynomial(&quotient_term)
            ));
            lines.push(format!(
                "- subtract: {}",
                format_dense_polynomial(&subtraction_term)
            ));
            lines.push(format!(
                "- new remainder: {}",
                format_dense_polynomial(&next_remainder)
            ));

            running_remainder = next_remainder;
            step += 1;
        }
    }

    lines.push(format!(
        "final quotient: {}",
        format_dense_polynomial(&quotient)
    ));
    lines.push(format!(
        "final remainder: {}",
        format_dense_polynomial(&remainder)
    ));

    let recomposed = divisor.mul(&quotient).add(&remainder);
    lines.push(format!(
        "verification: {} * {} + {} = {}",
        format_dense_polynomial(divisor),
        format_dense_polynomial(&quotient),
        format_dense_polynomial(&remainder),
        format_dense_polynomial(&recomposed)
    ));

    Ok(lines.join("\n"))
}

/// Builds the monomial polynomial `coefficient * x^degree`.
fn build_monomial<F: Field>(coefficient: F::Elem, degree: usize) -> DensePolynomial<F> {
    if F::is_zero(&coefficient) {
        return DensePolynomial::new(Vec::new());
    }

    let mut coefficients = vec![F::zero(); degree];
    coefficients.push(coefficient);
    DensePolynomial::new(coefficients)
}

/// Shifts a dense polynomial by multiplying it by `x^degree`.
fn shift_dense<F: Field>(polynomial: &DensePolynomial<F>, degree: usize) -> DensePolynomial<F> {
    if polynomial.is_zero() {
        return DensePolynomial::new(Vec::new());
    }

    let mut coefficients = vec![F::zero(); degree];
    coefficients.extend(polynomial.coefficients().iter().cloned());
    DensePolynomial::new(coefficients)
}

#[cfg(test)]
mod tests {
    use crate::fields::{Field, Fp};
    use crate::polynomials::{DensePolynomial, PolynomialError};

    use super::explain_dense_division;

    type F17 = Fp<17>;

    fn dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    #[test]
    fn dense_division_explanation_shows_exact_division_steps() {
        let explanation = explain_dense_division(&dense(&[2, 3, 1]), &dense(&[1, 1]))
            .expect("division should work");

        assert!(explanation.contains("Dense polynomial division"));
        assert!(explanation.contains("step 0:"));
        assert!(explanation.contains("step 1:"));
        assert!(explanation.contains("final quotient: x + 2 (mod 17)"));
        assert!(explanation.contains("final remainder: 0"));
    }

    #[test]
    fn dense_division_explanation_shows_non_zero_remainder() {
        let explanation = explain_dense_division(&dense(&[1, 2, 0, 1]), &dense(&[1, 0, 1]))
            .expect("division should work");

        assert!(explanation.contains("final quotient: x"));
        assert!(explanation.contains("final remainder: x + 1 (mod 17)"));
        assert!(explanation.contains("verification: x^2 + 1 (mod 17) * x + x + 1 (mod 17)"));
    }

    #[test]
    fn dense_division_explanation_handles_high_degree_divisor() {
        let explanation = explain_dense_division(&dense(&[3, 5]), &dense(&[1, 0, 1]))
            .expect("division should work");

        assert!(explanation.contains("divisor degree is larger than the dividend degree"));
        assert!(explanation.contains("final quotient: 0"));
        assert!(explanation.contains("final remainder: 5 (mod 17)*x + 3 (mod 17)"));
    }

    #[test]
    fn dense_division_explanation_rejects_zero_divisor() {
        let error =
            explain_dense_division(&dense(&[1, 2, 3]), &DensePolynomial::<F17>::new(Vec::new()))
                .expect_err("zero divisor should fail");

        assert_eq!(error, PolynomialError::DivisionByZeroPolynomial);
    }
}
