use crate::fields::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};
use crate::visualization::VisualizableField;

use super::format_dense_polynomial;

/// Explains the Euclidean algorithm for dense univariate polynomials over a
/// field.
///
/// The explanation shows the successive remainder steps
///
/// `gcd(a, b) = gcd(b, a mod b)`
///
/// and finishes by reporting the monic gcd chosen by the library.
pub fn explain_dense_gcd<F>(
    lhs: &DensePolynomial<F>,
    rhs: &DensePolynomial<F>,
) -> Result<String, PolynomialError>
where
    F: Field,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Dense polynomial gcd".to_string(),
        format!("lhs: {}", format_dense_polynomial(lhs)),
        format!("rhs: {}", format_dense_polynomial(rhs)),
        "method: Euclidean algorithm on successive remainders".to_string(),
    ];

    if lhs.is_zero() && rhs.is_zero() {
        lines.push(
            "both inputs are zero, so the gcd is defined here as the zero polynomial".to_string(),
        );
        lines.push("result: 0".to_string());
        return Ok(lines.join("\n"));
    }

    let mut a = lhs.clone();
    let mut b = rhs.clone();
    let mut step = 0usize;

    while !a.is_zero() && !b.is_zero() {
        let remainder = a.rem(&b)?;
        lines.push(format!(
            "step {step}: {} mod {} = {}",
            format_dense_polynomial(&a),
            format_dense_polynomial(&b),
            format_dense_polynomial(&remainder)
        ));
        a = b;
        b = remainder;
        step += 1;
    }

    let gcd = lhs.gcd(rhs);
    lines.push(format!("monic gcd: {}", format_dense_polynomial(&gcd)));
    Ok(lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use crate::fields::{Field, Fp};
    use crate::polynomials::DensePolynomial;

    use super::explain_dense_gcd;

    type F17 = Fp<17>;

    fn dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    #[test]
    fn dense_gcd_explanation_shows_remainder_chain_and_result() {
        let lhs = dense(&[2, 3, 1]);
        let rhs = dense(&[1, 3, 3, 1]);
        let explanation = explain_dense_gcd(&lhs, &rhs).expect("gcd should work");

        assert!(explanation.contains("Dense polynomial gcd"));
        assert!(explanation.contains("step 0:"));
        assert!(explanation.contains("monic gcd: x + 1 (mod 17)"));
    }

    #[test]
    fn dense_gcd_explanation_handles_double_zero() {
        let zero = DensePolynomial::<F17>::new(Vec::new());
        let explanation = explain_dense_gcd(&zero, &zero).expect("gcd should work");

        assert!(explanation.contains("both inputs are zero"));
        assert!(explanation.contains("result: 0"));
    }
}
