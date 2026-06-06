use crate::fields::Field;
use crate::polynomials::{
    DensePolynomial, IrreducibilityBackend, IrreducibilityStatus, PolynomialError,
    ReducibilityReason, irreducibility_status,
};
use crate::visualization::VisualizableField;

use crate::visualization::polynomials::format_dense_polynomial;

/// Returns a short educational description of an irreducibility status.
pub fn describe_irreducibility_status<F>(status: &IrreducibilityStatus<F>) -> String
where
    F: Field,
    F::Elem: VisualizableField,
{
    match status {
        IrreducibilityStatus::Constant => {
            "constant polynomial; not irreducible by the library's convention".to_string()
        }
        IrreducibilityStatus::Linear => {
            "linear polynomial; irreducible by the library's convention".to_string()
        }
        IrreducibilityStatus::Irreducible => {
            "irreducible; no non-trivial factor was found by the current backend".to_string()
        }
        IrreducibilityStatus::Reducible { divisor, quotient } => format!(
            "reducible; witness factorization: {} = ({}) * ({})",
            format_dense_polynomial(&divisor.mul(quotient)),
            format_dense_polynomial(divisor),
            format_dense_polynomial(quotient)
        ),
        IrreducibilityStatus::ReducibleWithoutWitness { reason } => format!(
            "reducible; current backend conclusion: {}",
            format_reducibility_reason(*reason)
        ),
    }
}

/// Explains irreducibility classification for dense univariate polynomials.
///
/// This helper keeps the algorithmic explanation in the `polynomials`
/// visualization layer, separate from the field-specific interpretation of
/// quotient moduli handled under `visualization/fields/`.
///
/// Unlike most other explainers, this function treats
/// [`PolynomialError::UndeterminedIrreducibility`] as an educational outcome
/// rather than a hard failure: it returns a textual explanation of the
/// inconclusive exact partial backend instead of bubbling the error up.
pub fn explain_dense_irreducibility<F>(
    polynomial: &DensePolynomial<F>,
) -> Result<String, PolynomialError>
where
    F: Field + IrreducibilityBackend,
    F::Elem: VisualizableField,
{
    let mut lines = vec![
        "Dense polynomial irreducibility".to_string(),
        format!("polynomial: {}", format_dense_polynomial(polynomial)),
        format!(
            "base field algebraically closed: {}",
            if F::IS_ALGEBRAICALLY_CLOSED {
                "yes"
            } else {
                "no"
            }
        ),
    ];

    match irreducibility_status(polynomial) {
        Ok(status) => {
            lines.push(format!(
                "status summary: {}",
                describe_irreducibility_status(&status)
            ));

            match status {
                IrreducibilityStatus::Constant => {
                    lines.push("convention: constants are not considered irreducible".to_string());
                }
                IrreducibilityStatus::Linear => {
                    lines.push(
                        "convention: linear polynomials are considered irreducible".to_string(),
                    );
                }
                IrreducibilityStatus::Irreducible => {
                    lines.push("conclusion: the current backend certifies that no non-trivial factor is required to explain the polynomial".to_string());
                }
                IrreducibilityStatus::Reducible { divisor, quotient } => {
                    lines.push(format!(
                        "witness divisor: {}",
                        format_dense_polynomial(&divisor)
                    ));
                    lines.push(format!(
                        "witness quotient: {}",
                        format_dense_polynomial(&quotient)
                    ));
                    lines.push(format!(
                        "factorization: {} = ({}) * ({})",
                        format_dense_polynomial(polynomial),
                        format_dense_polynomial(&divisor),
                        format_dense_polynomial(&quotient)
                    ));
                }
                IrreducibilityStatus::ReducibleWithoutWitness { reason } => {
                    lines.push(format!("reason: {}", format_reducibility_reason(reason)));
                    lines.push(
                        "note: the backend concluded reducibility without exposing a concrete factorization witness".to_string(),
                    );
                }
            }
        }
        Err(PolynomialError::UndeterminedIrreducibility(message)) => {
            lines.push("status summary: inconclusive".to_string());
            lines.push(
                "conclusion: the current backend is exact but partial for this field family"
                    .to_string(),
            );
            lines.push(format!("current backend message: {message}"));
        }
        Err(error) => return Err(error),
    }

    Ok(lines.join("\n"))
}

fn format_reducibility_reason(reason: ReducibilityReason) -> &'static str {
    match reason {
        ReducibilityReason::AlgebraicallyClosedField => {
            "the base field is algebraically closed, so every degree >= 2 polynomial factors non-trivially"
        }
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use crate::fields::{ComplexApprox, Field, Fp, Q};
    use crate::polynomials::{DensePolynomial, IrreducibilityStatus, ReducibilityReason};

    use crate::visualization::polynomials::{describe_irreducibility_status, explain_dense_irreducibility};

    type F17 = Fp<17>;

    fn dense(values: &[u64]) -> DensePolynomial<F17> {
        DensePolynomial::new(values.iter().copied().map(F17::elem_from_u64).collect())
    }

    fn complex_dense(values: &[(f64, f64)]) -> DensePolynomial<ComplexApprox> {
        DensePolynomial::new(
            values
                .iter()
                .copied()
                .map(|(re, im)| Complex64::new(re, im))
                .collect(),
        )
    }

    #[test]
    fn irreducibility_status_description_handles_witnessless_reducibility() {
        let status = IrreducibilityStatus::<ComplexApprox>::ReducibleWithoutWitness {
            reason: ReducibilityReason::AlgebraicallyClosedField,
        };

        let description = describe_irreducibility_status(&status);
        assert!(description.contains("reducible"));
        assert!(description.contains("algebraically closed"));
    }

    #[test]
    fn dense_irreducibility_explanation_reports_factorization_witness() {
        let explanation =
            explain_dense_irreducibility(&dense(&[1, 0, 1])).expect("explanation should work");

        assert!(explanation.contains("Dense polynomial irreducibility"));
        assert!(explanation.contains("status summary: reducible"));
        assert!(explanation.contains("witness divisor:"));
        assert!(explanation.contains("factorization:"));
    }

    #[test]
    fn dense_irreducibility_explanation_reports_irreducible_case() {
        let explanation =
            explain_dense_irreducibility(&dense(&[3, 0, 1])).expect("explanation should work");

        assert!(explanation.contains("status summary: irreducible"));
        assert!(explanation.contains("current backend certifies"));
    }

    #[test]
    fn dense_irreducibility_explanation_reports_theoretical_complex_reducibility() {
        let explanation =
            explain_dense_irreducibility(&complex_dense(&[(1.0, 0.0), (0.0, 0.0), (1.0, 0.0)]))
                .expect("explanation should work");

        assert!(explanation.contains("base field algebraically closed: yes"));
        assert!(explanation.contains("status summary: reducible"));
        assert!(explanation.contains("note: the backend concluded reducibility"));
    }

    #[test]
    fn dense_irreducibility_explanation_reports_inconclusive_q_cases() {
        use num_bigint::BigInt;
        use num_rational::BigRational;
        use num_traits::One;

        let leading = [2_u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
            .into_iter()
            .fold(BigInt::one(), |accumulator, prime| {
                accumulator * BigInt::from(prime)
            });
        let constant = &leading + BigInt::one();
        let polynomial = DensePolynomial::<Q>::new(vec![
            BigRational::from_integer(constant),
            Q::zero(),
            Q::zero(),
            Q::zero(),
            BigRational::from_integer(leading),
        ]);

        let explanation = explain_dense_irreducibility(&polynomial)
            .expect("inconclusive explanation should work");

        assert!(explanation.contains("status summary: inconclusive"));
        assert!(explanation.contains("exact but partial"));
        assert!(explanation.contains("rational-root search"));
    }
}
