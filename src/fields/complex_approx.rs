use crate::fields::traits::*;
use num_bigint::BigInt;
use num_complex::Complex64;
use num_traits::ToPrimitive;

use crate::fields::{
    FieldCharacteristic,
    error::FieldError,
    traits::{CbrtField, SqrtField},
};
use crate::numerics::ApproxTolerance;

/// Structured report for one approximate complex comparison.
///
/// Even though this report belongs to the [`ComplexApprox`] backend, the
/// stored values are concrete [`Complex64`] elements. `ComplexApprox` itself is
/// the field family marker, not the runtime value representation.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ApproxComparisonReport {
    /// Left-hand complex value that was compared.
    pub lhs: Complex64,
    /// Right-hand complex value that was compared.
    pub rhs: Complex64,
    /// Difference `lhs - rhs`.
    pub difference: Complex64,
    /// Euclidean norm of `lhs - rhs`.
    pub absolute_error: f64,
    /// Tolerance policy used for the comparison.
    pub tolerance: ApproxTolerance,
    /// Whether the values were considered close under `tolerance`.
    pub is_close: bool,
}

/// Approximate complex-number field backed by [`Complex64`].
///
/// This implementation is intended for numerical experimentation and testing,
/// not for exact algebraic reasoning. Equality and zero checks are approximate
/// and use [`ComplexApprox::default_tolerance()`] by default.
#[derive(Clone, Copy, Debug)]
pub struct ComplexApprox;

impl ComplexApprox {
    /// Returns the default tolerance policy for this approximate field model.
    ///
    /// The current default matches the crate's stricter educational preset so
    /// existing comparisons stay conservative while the tolerance policy moves
    /// out of a hardcoded scalar constant and into a reusable shared type.
    pub fn default_tolerance() -> ApproxTolerance {
        ApproxTolerance::strict()
    }

    /// Compares two complex values using an explicit tolerance policy.
    /// The comparison uses a mixed absolute/relative rule:
    /// `|x - y| <= max(abs_tol, rel_tol * max(|x|, |y|))`
    ///
    /// This keeps small values stable near zero while still allowing a
    /// scale-aware comparison for larger magnitudes.
    pub fn eq_with_tolerance(x: &Complex64, y: &Complex64, tolerance: ApproxTolerance) -> bool {
        Self::comparison_report(x, y, tolerance).is_close
    }

    /// Returns whether `x` is approximately zero under an explicit tolerance.
    ///
    /// Zero checks use only the absolute tolerance because there is no larger
    /// reference scale to compare against meaningfully.
    pub fn is_zero_with_tolerance(x: &Complex64, tolerance: ApproxTolerance) -> bool {
        x.norm() <= tolerance.absolute
    }

    /// Builds a structured comparison report under an explicit tolerance.
    ///
    /// The report records both compared values, their difference, the absolute
    /// error `|lhs - rhs|`, the tolerance policy used, and the final closeness
    /// verdict under the mixed absolute/relative comparison rule.
    pub(crate) fn comparison_report(
        lhs: &Complex64,
        rhs: &Complex64,
        tolerance: ApproxTolerance,
    ) -> ApproxComparisonReport {
        let difference = *lhs - *rhs;
        let absolute_error = difference.norm();
        let scale = lhs.norm().max(rhs.norm());
        let bound = tolerance.absolute.max(tolerance.relative * scale);

        ApproxComparisonReport {
            lhs: *lhs,
            rhs: *rhs,
            difference,
            absolute_error,
            tolerance,
            is_close: absolute_error <= bound,
        }
    }
}

impl Field for ComplexApprox {
    /// The complex numbers are algebraically closed.
    ///
    /// `ComplexApprox` is only an approximate numerical model of `C`, but it
    /// is still meant to advertise the mathematical closedness property of the
    /// backend field family it is approximating.
    const IS_ALGEBRAICALLY_CLOSED: bool = true;

    type Elem = Complex64;

    fn characteristic() -> FieldCharacteristic {
        FieldCharacteristic::Zero
    }

    /// Returns the additive identity.
    fn zero() -> Self::Elem {
        Self::Elem::new(0.0, 0.0)
    }

    /// Returns the multiplicative identity.
    fn one() -> Self::Elem {
        Self::Elem::new(1.0, 0.0)
    }

    /// Embeds an integer into the complex numbers, approximately.
    fn from_bigint(n: &BigInt) -> Self::Elem {
        Self::Elem::new(
            n.to_f64()
                .expect("ComplexApprox integer literals should fit in f64 for this route"),
            0.0,
        )
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
        if Self::is_zero_with_tolerance(x, Self::default_tolerance()) {
            None
        } else {
            Some(Self::one() / *x)
        }
    }

    /// Returns whether two complex values are equal under the default tolerance.
    fn eq(x: &Self::Elem, y: &Self::Elem) -> bool {
        Self::eq_with_tolerance(x, y, Self::default_tolerance())
    }

    /// Returns the multiplicative inverse or a structured zero-division error.
    fn inverse(x: &Self::Elem) -> Result<Self::Elem, FieldError> {
        Self::inv(x).ok_or(FieldError::DivisionByZero)
    }
}

impl SqrtField for ComplexApprox {
    /// Returns the principal complex square root from the numerical backend.
    ///
    /// This is an approximate floating-point computation, not an exact algebraic
    /// square-root procedure.
    fn sqrt(x: &Self::Elem) -> Option<Self::Elem> {
        Some(x.sqrt())
    }
}

impl CbrtField for ComplexApprox {
    /// Returns the principal complex cube root from the numerical backend.
    ///
    /// This is an approximate floating-point computation, not an exact
    /// algebraic cube-root procedure. As with the principal complex square
    /// root, this value is branch-sensitive.
    fn cbrt(x: &Self::Elem) -> Option<Self::Elem> {
        Some(x.powf(1.0 / 3.0))
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::traits::*;
    use std::hint::black_box;

    use num_bigint::BigUint;
    use num_complex::Complex64;

    use crate::{
        fields::{
            ComplexApprox, FieldError,
            complex_approx::ApproxComparisonReport,
            traits::{CbrtField, SqrtField},
        },
        numerics::ApproxTolerance,
    };

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
        assert!(ComplexApprox::has_characteristic(0));
        assert_close(ComplexApprox::zero(), c(0.0, 0.0));
        assert_close(ComplexApprox::one(), c(1.0, 0.0));
        assert_close(ComplexApprox::from_i64(-7), c(-7.0, 0.0));
        assert_close(ComplexApprox::from_i64(42), c(42.0, 0.0));
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
        assert_close(
            ComplexApprox::pow(&i, &BigUint::from(0u8)),
            ComplexApprox::one(),
        );
        assert_close(
            ComplexApprox::pow(&i, &BigUint::from(4u8)),
            ComplexApprox::one(),
        );
        assert_close(ComplexApprox::pow(&z, &BigUint::from(3u8)), z * z * z);
    }

    #[test]
    fn approximate_equality_uses_the_default_tolerance() {
        let tolerance = ComplexApprox::default_tolerance();
        let x = c(1.0, -2.0);
        let inside = c(1.0 + tolerance.absolute * 0.25, -2.0);
        let outside = c(1.0 + tolerance.absolute * 4.0, -2.0);

        assert!(ComplexApprox::eq(&x, &inside));
        assert!(!ComplexApprox::eq(&x, &outside));
    }

    #[test]
    fn zero_check_uses_the_default_absolute_tolerance() {
        let tolerance = ComplexApprox::default_tolerance();
        let tiny = c(tolerance.absolute * 0.2, -tolerance.absolute * 0.2);
        let not_tiny = c(tolerance.absolute * 10.0, 0.0);

        assert!(ComplexApprox::is_zero(&tiny));
        assert!(!ComplexApprox::is_zero(&not_tiny));
    }

    #[test]
    fn explicit_tolerance_can_relax_equality() {
        let x = c(1000.0, 0.0);
        let y = c(1000.0005, 0.0);
        let loose = ApproxTolerance::new(1.0e-9, 1.0e-6);

        assert!(ComplexApprox::eq_with_tolerance(&x, &y, loose));
        assert!(!ComplexApprox::eq(&x, &y));
    }

    #[test]
    fn comparison_report_records_operands_difference_and_verdict() {
        let lhs = c(3.0, 4.0);
        let rhs = c(3.0, 4.0000000000005);
        let tolerance = ApproxTolerance::strict();
        let report = ComplexApprox::comparison_report(&lhs, &rhs, tolerance);

        assert_eq!(
            report,
            ApproxComparisonReport {
                lhs,
                rhs,
                difference: lhs - rhs,
                absolute_error: (lhs - rhs).norm(),
                tolerance,
                is_close: true,
            }
        );
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

    #[test]
    fn algebraic_closedness_metadata_matches_complex_numbers() {
        assert!(black_box(ComplexApprox::IS_ALGEBRAICALLY_CLOSED));
    }

    #[test]
    fn sqrt_returns_a_complex_root_for_negative_reals() {
        let root = ComplexApprox::sqrt(&c(-1.0, 0.0)).expect("complex numbers have square roots");

        assert_close(ComplexApprox::square(&root), c(-1.0, 0.0));
    }

    #[test]
    fn sqrt_pair_returns_opposite_complex_roots() {
        let (left, right) =
            ComplexApprox::sqrt_pair(&c(4.0, 0.0)).expect("complex numbers have square roots");

        assert_close(ComplexApprox::square(&left), c(4.0, 0.0));
        assert_close(ComplexApprox::square(&right), c(4.0, 0.0));
        assert_close(right, ComplexApprox::neg(&left));
    }

    #[test]
    fn cbrt_returns_principal_branch_value() {
        let root = ComplexApprox::cbrt(&c(-1.0, 0.0)).expect("complex numbers have cube roots");

        assert_close(ComplexApprox::cube(&root), c(-1.0, 0.0));
        assert!(root.im > 0.0);
    }

    #[test]
    fn default_tolerance_matches_the_shared_strict_preset() {
        assert_eq!(
            ComplexApprox::default_tolerance(),
            ApproxTolerance::strict()
        );
    }
}
