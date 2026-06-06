use num_complex::Complex64;
use num_rational::BigRational;

use crate::fields::{
    ComplexApprox, errors::FieldError, prime_field::Fp, rationals::Q, sqrt_field::SqrtField,
    traits::Field,
};
use crate::visualization::fields::{
    complex_approx::format_complex, prime_field::format_fp_elem, rationals::format_rational,
};

/// Explains how the current prime-field square-root backend behaves.
///
/// The explanation is intentionally honest about scope:
///
/// - `GF(2)` is handled as a tiny special case
/// - odd primes use Tonelli-Shanks
/// - the report says explicitly whether the input is a quadratic residue
pub fn explain_prime_field_square_root<const P: u64>(value: u64) -> Result<String, FieldError> {
    Fp::<P>::validate_modulus()?;
    let element = Fp::<P>::new_elem(value)?;
    let root_pair = Fp::<P>::sqrt_pair(&element);

    let algorithm = if P == 2 {
        "special case for GF(2)"
    } else {
        "Tonelli-Shanks over an odd prime field"
    };

    let mut lines = vec![
        format!("Square roots in GF({P})"),
        format!("input: {}", format_fp_elem(&element)),
        format!("algorithm: {algorithm}"),
    ];

    match root_pair {
        Some((left, right)) => {
            let left_check = Fp::<P>::square(&left);
            let right_check = Fp::<P>::square(&right);
            lines.push("quadratic residue: yes".to_string());
            lines.push(format!(
                "root pair: {}, {}",
                format_fp_elem(&left),
                format_fp_elem(&right)
            ));
            lines.push(format!(
                "verification: {}^2 = {}, {}^2 = {}",
                format_fp_elem(&left),
                format_fp_elem(&left_check),
                format_fp_elem(&right),
                format_fp_elem(&right_check)
            ));
        }
        None => {
            lines.push("quadratic residue: no".to_string());
            lines.push("result: no square root exists in this prime field".to_string());
        }
    }

    Ok(lines.join("\n"))
}

/// Explains rational square roots with exact arithmetic.
///
/// The current backend succeeds only when the reduced numerator and
/// denominator are both perfect integer squares.
pub fn explain_rational_square_root(x: &BigRational) -> String {
    let mut lines = vec![
        "Square roots in Q".to_string(),
        format!("input: {}", format_rational(x)),
        "current exact scope: succeeds only for rational squares already present in Q".to_string(),
    ];

    match Q::sqrt_pair(x) {
        Some((left, right)) => {
            lines.push("square root exists in Q: yes".to_string());
            lines.push(format!(
                "root pair: {}, {}",
                format_rational(&left),
                format_rational(&right)
            ));
            lines.push(format!(
                "verification: ({})^2 = {}",
                format_rational(&left),
                format_rational(&Q::square(&left))
            ));
        }
        None => {
            lines.push("square root exists in Q: no".to_string());
            lines.push("result: no exact rational square root".to_string());
        }
    }

    lines.join("\n")
}

/// Explains approximate complex square roots.
///
/// `ComplexApprox` returns the principal square root from the numerical
/// complex backend and obtains the other root by negation.
pub fn explain_complex_square_root(z: &Complex64) -> String {
    let (principal, opposite) =
        ComplexApprox::sqrt_pair(z).expect("complex numbers always admit square roots");

    format!(
        "Square roots in C (approx)\n\
         input: {}\n\
         branch: principal square root from the numerical backend\n\
         root pair: {}, {}\n\
         verification: ({})^2 = {}",
        format_complex(z),
        format_complex(&principal),
        format_complex(&opposite),
        format_complex(&principal),
        format_complex(&ComplexApprox::square(&principal))
    )
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_complex::Complex64;
    use num_rational::BigRational;

    use crate::visualization::fields::{
        explain_complex_square_root, explain_prime_field_square_root, explain_rational_square_root,
    };

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn prime_field_square_root_explanation_shows_residue_case() {
        let explanation =
            explain_prime_field_square_root::<17>(4).expect("prime field explanation should work");

        assert!(explanation.contains("Square roots in GF(17)"));
        assert!(explanation.contains("algorithm: Tonelli-Shanks"));
        assert!(explanation.contains("quadratic residue: yes"));
        assert!(explanation.contains("input: 4"));
    }

    #[test]
    fn prime_field_square_root_explanation_shows_non_residue_case() {
        let explanation =
            explain_prime_field_square_root::<17>(3).expect("prime field explanation should work");

        assert!(explanation.contains("quadratic residue: no"));
        assert!(explanation.contains("no square root exists"));
    }

    #[test]
    fn rational_square_root_explanation_is_honest_about_exact_scope() {
        let explanation = explain_rational_square_root(&q(9, 16));

        assert!(explanation.contains("Square roots in Q"));
        assert!(explanation.contains("current exact scope"));
        assert!(explanation.contains("root pair: 3/4, -3/4"));
        assert!(explanation.contains("verification: (3/4)^2 = 9/16"));
    }

    #[test]
    fn rational_square_root_explanation_reports_absence_of_exact_root() {
        let explanation = explain_rational_square_root(&q(2, 1));

        assert!(explanation.contains("square root exists in Q: no"));
        assert!(explanation.contains("no exact rational square root"));
    }

    #[test]
    fn complex_square_root_explanation_shows_principal_branch() {
        let explanation = explain_complex_square_root(&Complex64::new(-1.0, 0.0));

        assert!(explanation.contains("Square roots in C (approx)"));
        assert!(explanation.contains("branch: principal square root"));
        assert!(explanation.contains("0.000000 + 1.000000i"));
    }
}
