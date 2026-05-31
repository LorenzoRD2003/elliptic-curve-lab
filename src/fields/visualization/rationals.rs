use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::fields::{
    errors::FieldError, rationals::Q, traits::Field, visualization::traits::Visualizable,
};

/// Returns a short textual description of the rational field `Q`.
pub fn format_rational_field() -> String {
    "Q\ncharacteristic: 0\nfinite: no\nexact arithmetic: yes".to_string()
}

/// Formats a rational number in canonical exact form.
///
/// Integers are shown without a denominator. Non-integral rationals are shown
/// as `n/d` with the sign carried by the numerator.
pub fn format_rational(x: &BigRational) -> String {
    if x.denom().is_one() {
        x.numer().to_string()
    } else {
        format!("{}/{}", x.numer(), x.denom())
    }
}

/// Returns a short educational description of a rational value.
pub fn describe_rational(x: &BigRational) -> String {
    let sign = if x.is_zero() {
        "zero"
    } else if x.is_positive() {
        "positive"
    } else {
        "negative"
    };

    format!(
        "q = {}\nnumerator: {}\ndenominator: {}\nsign: {}\nintegral: {}\nzero: {}",
        format_rational(x),
        x.numer(),
        x.denom(),
        sign,
        x.denom().is_one(),
        x.is_zero()
    )
}

/// Explains an exact addition in `Q`.
pub fn explain_rational_add(lhs: &BigRational, rhs: &BigRational) -> String {
    let a = lhs.numer();
    let b = lhs.denom();
    let c = rhs.numer();
    let d = rhs.denom();
    let left_scaled = a * d;
    let right_scaled = c * b;
    let denominator = b * d;
    let result = Q::add(lhs, rhs);

    format!(
        "Addition in Q\n\
         lhs: {}\n\
         rhs: {}\n\
         common denominator: {} * {} = {}\n\
         scaled numerators: ({} * {}) + ({} * {}) = {} + {} = {}\n\
         result: {}",
        format_rational(lhs),
        format_rational(rhs),
        b,
        d,
        denominator,
        a,
        d,
        c,
        b,
        left_scaled,
        right_scaled,
        left_scaled.clone() + right_scaled.clone(),
        format_rational(&result)
    )
}

/// Explains an exact multiplication in `Q`.
pub fn explain_rational_mul(lhs: &BigRational, rhs: &BigRational) -> String {
    let a = lhs.numer();
    let b = lhs.denom();
    let c = rhs.numer();
    let d = rhs.denom();
    let numerator = a * c;
    let denominator = b * d;
    let result = Q::mul(lhs, rhs);

    format!(
        "Multiplication in Q\n\
         lhs: {}\n\
         rhs: {}\n\
         numerator: {} * {} = {}\n\
         denominator: {} * {} = {}\n\
         result: {}",
        format_rational(lhs),
        format_rational(rhs),
        a,
        c,
        numerator,
        b,
        d,
        denominator,
        format_rational(&result)
    )
}

/// Explains the multiplicative inverse of a non-zero rational.
pub fn explain_rational_inverse(x: &BigRational) -> Result<String, FieldError> {
    let inverse = Q::inverse(x)?;

    Ok(format!(
        "Inverse in Q\n\
         element: {}\n\
         inverse: {}\n\
         verification: {} * {} = {}",
        format_rational(x),
        format_rational(&inverse),
        format_rational(x),
        format_rational(&inverse),
        format_rational(&Q::mul(x, &inverse))
    ))
}

/// Explains exact division in `Q`.
pub fn explain_rational_div(lhs: &BigRational, rhs: &BigRational) -> Result<String, FieldError> {
    let reciprocal = Q::inverse(rhs)?;
    let result = Q::div(lhs, rhs)?;

    Ok(format!(
        "Division in Q\n\
         lhs: {}\n\
         rhs: {}\n\
         reciprocal of rhs: {}\n\
         reduction to multiplication: {} * {} = {}",
        format_rational(lhs),
        format_rational(rhs),
        format_rational(&reciprocal),
        format_rational(lhs),
        format_rational(&reciprocal),
        format_rational(&result)
    ))
}

impl Visualizable for BigRational {
    fn format_elem(&self) -> String {
        format_rational(self)
    }

    fn describe(&self) -> String {
        describe_rational(self)
    }

    fn inverse(&self) -> Option<String> {
        Q::inv(self).map(|value| format_rational(&value))
    }

    fn explain_add(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_rational_add(lhs, rhs))
    }

    fn explain_mul(lhs: &Self, rhs: &Self) -> Option<String> {
        Some(explain_rational_mul(lhs, rhs))
    }

    fn explain_div(lhs: &Self, rhs: &Self) -> Option<String> {
        explain_rational_div(lhs, rhs).ok()
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;
    use num_rational::BigRational;

    use super::{
        describe_rational, explain_rational_add, explain_rational_div, explain_rational_inverse,
        explain_rational_mul, format_rational, format_rational_field,
    };
    use crate::fields::FieldError;
    use crate::fields::visualization::Visualizable;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn rational_field_summary_is_readable() {
        let summary = format_rational_field();
        assert!(summary.contains("Q"));
        assert!(summary.contains("characteristic: 0"));
        assert!(summary.contains("finite: no"));
        assert!(summary.contains("exact arithmetic: yes"));
    }

    #[test]
    fn rational_formatting_is_canonical_and_compact() {
        assert_eq!(format_rational(&q(6, 3)), "2");
        assert_eq!(format_rational(&q(-4, 6)), "-2/3");
        assert_eq!(format_rational(&q(0, 5)), "0");
    }

    #[test]
    fn rational_description_contains_key_metadata() {
        let description = describe_rational(&q(-7, 3));
        assert!(description.contains("q = -7/3"));
        assert!(description.contains("numerator: -7"));
        assert!(description.contains("denominator: 3"));
        assert!(description.contains("sign: negative"));
        assert!(description.contains("integral: false"));
    }

    #[test]
    fn rational_addition_explanation_shows_common_denominator() {
        let explanation = explain_rational_add(&q(2, 3), &q(5, 7));
        assert!(explanation.contains("lhs: 2/3"));
        assert!(explanation.contains("rhs: 5/7"));
        assert!(explanation.contains("common denominator: 3 * 7 = 21"));
        assert!(explanation.contains("scaled numerators: (2 * 7) + (5 * 3) = 14 + 15 = 29"));
        assert!(explanation.contains("result: 29/21"));
    }

    #[test]
    fn rational_multiplication_explanation_is_exact() {
        let explanation = explain_rational_mul(&q(2, 3), &q(5, 7));
        assert!(explanation.contains("numerator: 2 * 5 = 10"));
        assert!(explanation.contains("denominator: 3 * 7 = 21"));
        assert!(explanation.contains("result: 10/21"));
    }

    #[test]
    fn rational_inverse_explanation_shows_verification() {
        let explanation = explain_rational_inverse(&q(-3, 5)).expect("inverse should exist");
        assert!(explanation.contains("element: -3/5"));
        assert!(explanation.contains("inverse: -5/3"));
        assert!(explanation.contains("verification: -3/5 * -5/3 = 1"));
    }

    #[test]
    fn rational_inverse_rejects_zero() {
        let error = explain_rational_inverse(&q(0, 1)).expect_err("zero has no inverse");
        assert!(matches!(error, FieldError::DivisionByZero));
    }

    #[test]
    fn rational_division_explanation_reduces_to_multiplication() {
        let explanation = explain_rational_div(&q(2, 3), &q(5, 7)).expect("division should work");
        assert!(explanation.contains("lhs: 2/3"));
        assert!(explanation.contains("rhs: 5/7"));
        assert!(explanation.contains("reciprocal of rhs: 7/5"));
        assert!(explanation.contains("reduction to multiplication: 2/3 * 7/5 = 14/15"));
    }

    #[test]
    fn rational_division_rejects_zero_denominator_operand() {
        let error = explain_rational_div(&q(1, 2), &q(0, 1)).expect_err("division by zero fails");
        assert!(matches!(error, FieldError::DivisionByZero));
    }

    #[test]
    fn rational_visualizable_trait_reuses_core_helpers() {
        let lhs = q(2, 3);
        let rhs = q(5, 7);

        assert_eq!(lhs.format_compact(), "2/3");
        assert_eq!(lhs.format_elem(), "2/3");
        assert!(lhs.describe().contains("numerator: 2"));
        assert_eq!(lhs.inverse().expect("inverse should exist"), "3/2");

        let add =
            BigRational::explain_add(&lhs, &rhs).expect("rational addition should be explainable");
        assert!(add.contains("Addition in Q"));
        assert!(add.contains("result: 29/21"));

        let mul = BigRational::explain_mul(&lhs, &rhs)
            .expect("rational multiplication should be explainable");
        assert!(mul.contains("Multiplication in Q"));
        assert!(mul.contains("result: 10/21"));

        let div =
            BigRational::explain_div(&lhs, &rhs).expect("rational division should be explainable");
        assert!(div.contains("Division in Q"));
        assert!(div.contains("14/15"));
    }
}
