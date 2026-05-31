use num_complex::Complex64;

use crate::fields::{errors::FieldError, traits::Field};

/// Approximate complex-number field backed by [`Complex64`].
///
/// This implementation is intended for numerical experimentation and testing,
/// not for exact algebraic reasoning. Equality and zero checks are approximate
/// and use [`ComplexApprox::EPS`] as a tolerance.
#[derive(Clone, Copy, Debug)]
pub struct ComplexApprox;

impl ComplexApprox {
    /// Absolute tolerance used for approximate equality and zero checks.
    pub const EPS: f64 = 1.0e-12;
}

impl Field for ComplexApprox {
    type Elem = Complex64;

    /// Returns the additive identity.
    fn zero() -> Self::Elem {
        Self::Elem::new(0.0, 0.0)
    }

    /// Returns the multiplicative identity.
    fn one() -> Self::Elem {
        Self::Elem::new(1.0, 0.0)
    }

    /// Embeds a signed integer into the complex numbers.
    fn from_i64(n: i64) -> Self::Elem {
        Self::Elem::new(n as f64, 0.0)
    }

    /// Adds two complex numbers.
    fn add(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        *x + *y
    }

    /// Subtracts `y` from `x`.
    fn sub(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        *x - *y
    }

    /// Multiplies two complex numbers.
    fn mul(x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        *x * *y
    }

    /// Returns the additive inverse of `x`.
    fn neg(x: &Self::Elem) -> Self::Elem {
        -*x
    }

    /// Returns the multiplicative inverse when `x` is not approximately zero.
    fn inv(x: &Self::Elem) -> Option<Self::Elem> {
        if Self::is_zero(x) {
            None
        } else {
            Some(Self::one() / *x)
        }
    }

    /// Returns whether two complex values are equal up to [`ComplexApprox::EPS`].
    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        (*x - *y).norm() <= Self::EPS
    }

    /// Returns the multiplicative inverse or a structured zero-division error.
    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::inv(x).ok_or(FieldError::DivisionByZero)
    }

    /// Embeds an unsigned integer into the complex numbers.
    fn elem_from_u64(value: u64) -> Self::Elem {
        Self::Elem::new(value as f64, 0.0)
    }
}

#[cfg(test)]
mod tests {
    use num_complex::Complex64;

    use super::ComplexApprox;
    use crate::fields::{Field, FieldError};

    fn c(re: f64, im: f64) -> Complex64 {
        Complex64::new(re, im)
    }

    fn assert_close(actual: Complex64, expected: Complex64) {
        assert!(
            ComplexApprox::eq(&actual, &expected),
            "expected {expected:?}, got {actual:?}"
        );
    }

    #[test]
    fn zero_one_and_integer_embeddings_are_correct() {
        assert_close(ComplexApprox::zero(), c(0.0, 0.0));
        assert_close(ComplexApprox::one(), c(1.0, 0.0));
        assert_close(ComplexApprox::from_i64(-7), c(-7.0, 0.0));
        assert_close(ComplexApprox::elem_from_u64(42), c(42.0, 0.0));
    }

    #[test]
    fn basic_arithmetic_works() {
        let x = c(2.0, 3.0);
        let y = c(-1.5, 0.5);

        assert_close(ComplexApprox::add(&x, &y), c(0.5, 3.5));
        assert_close(ComplexApprox::sub(&x, &y), c(3.5, 2.5));
        assert_close(ComplexApprox::mul(&x, &y), c(-4.5, -3.5));
        assert_close(ComplexApprox::neg(&x), c(-2.0, -3.0));
    }

    #[test]
    fn square_cube_and_pow_work() {
        let i = c(0.0, 1.0);
        let z = c(1.0, 2.0);

        assert_close(ComplexApprox::square(&i), c(-1.0, 0.0));
        assert_close(ComplexApprox::cube(&i), c(0.0, -1.0));
        assert_close(ComplexApprox::pow(&i, 0), ComplexApprox::one());
        assert_close(ComplexApprox::pow(&i, 4), ComplexApprox::one());
        assert_close(ComplexApprox::pow(&z, 3), z * z * z);
    }

    #[test]
    fn approximate_equality_uses_eps() {
        let x = c(1.0, -2.0);
        let inside = c(1.0 + ComplexApprox::EPS * 0.25, -2.0);
        let outside = c(1.0 + ComplexApprox::EPS * 4.0, -2.0);

        assert!(ComplexApprox::eq(&x, &inside));
        assert!(!ComplexApprox::eq(&x, &outside));
    }

    #[test]
    fn zero_check_is_approximate() {
        let tiny = c(ComplexApprox::EPS * 0.2, -ComplexApprox::EPS * 0.2);
        let not_tiny = c(ComplexApprox::EPS * 10.0, 0.0);

        assert!(ComplexApprox::is_zero(&tiny));
        assert!(!ComplexApprox::is_zero(&not_tiny));
    }

    #[test]
    fn inverse_behaves_as_expected() {
        let z = c(2.0, -1.0);
        let inverse = ComplexApprox::inv(&z).expect("non-zero value should be invertible");

        assert_close(ComplexApprox::mul(&z, &inverse), ComplexApprox::one());
        assert_close(ComplexApprox::mul(&inverse, &z), ComplexApprox::one());
    }

    #[test]
    fn inverse_and_division_report_zero_division() {
        let zero = ComplexApprox::zero();
        let x = c(3.0, 4.0);

        assert_eq!(ComplexApprox::inv(&zero), None);
        assert!(matches!(
            ComplexApprox::inverse(&zero),
            Err(FieldError::DivisionByZero)
        ));
        assert!(matches!(
            ComplexApprox::div(&x, &zero),
            Err(FieldError::DivisionByZero)
        ));
    }

    #[test]
    fn division_matches_manual_complex_division() {
        let x = c(3.0, 4.0);
        let y = c(1.0, -2.0);
        let quotient = ComplexApprox::div(&x, &y).expect("division should succeed");

        assert_close(quotient, x / y);
        assert_close(ComplexApprox::mul(&quotient, &y), x);
    }

    #[test]
    fn additive_and_multiplicative_identities_hold_on_samples() {
        let samples = [c(0.0, 0.0), c(1.5, -2.5), c(-3.0, 4.25)];

        for sample in samples {
            assert_close(ComplexApprox::add(&sample, &ComplexApprox::zero()), sample);
            assert_close(ComplexApprox::mul(&sample, &ComplexApprox::one()), sample);
        }
    }
}
