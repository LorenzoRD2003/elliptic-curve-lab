//! Shared floating-point comparison tolerances for educational numerical code.
//!
//! This module is intentionally sibling infrastructure: it can be used by
//! approximate field backends, complex-analytic elliptic-curve helpers, or any
//! other crate surface that needs small explicit floating-point policy objects.

use num_complex::Complex64;

/// Absolute and relative tolerances for approximate floating-point comparison.
///
/// The current project uses this type for educational numerical experiments,
/// not certified numerical analysis. It therefore keeps the policy surface
/// intentionally small and easy to explain.
///
/// `ApproxTolerance::new(...)` does not validate its arguments. It simply
/// stores caller-supplied tolerances for later algorithms, which may decide to
/// reject invalid or nonsensical settings explicitly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ApproxTolerance {
    /// Absolute error budget for near-zero comparisons.
    pub absolute: f64,
    /// Relative error budget for scale-sensitive comparisons.
    pub relative: f64,
}

impl ApproxTolerance {
    /// Builds a tolerance package from explicit absolute and relative bounds.
    pub fn new(absolute: f64, relative: f64) -> Self {
        Self { absolute, relative }
    }

    /// Compares two real values with this mixed absolute/relative tolerance.
    ///
    /// The rule is
    ///
    /// `|lhs - rhs| <= max(absolute, relative * max(|lhs|, |rhs|))`.
    ///
    /// This keeps near-zero comparisons stable through the absolute tolerance
    /// while still allowing scale-aware comparisons for larger magnitudes.
    pub fn real_close(self, lhs: f64, rhs: f64) -> bool {
        let scale = lhs.abs().max(rhs.abs());
        let bound = self.absolute.max(self.relative * scale);

        (lhs - rhs).abs() <= bound
    }

    /// Returns the baseline tolerance preset for educational experiments.
    pub fn educational_default() -> Self {
        Self::new(1.0e-9, 1.0e-9)
    }

    /// Returns a tighter preset for more delicate numerical comparisons.
    pub fn strict() -> Self {
        Self::new(1.0e-12, 1.0e-12)
    }

    /// Returns a more forgiving preset for coarse exploratory computations.
    pub fn loose() -> Self {
        Self::new(1.0e-6, 1.0e-6)
    }
}

pub(crate) fn infinity_proximity_scale(tolerance: ApproxTolerance) -> f64 {
    tolerance
        .absolute
        .max(tolerance.relative)
        .max(f64::EPSILON.sqrt())
}

pub(crate) fn reciprocal_singularity_threshold(tolerance: ApproxTolerance) -> f64 {
    1.0 / infinity_proximity_scale(tolerance)
}

pub(crate) fn projective_unit_singularity_distance(value: &Complex64) -> f64 {
    let norm = value.norm();
    if norm == 0.0 {
        return 0.0;
    }

    norm.min((Complex64::new(1.0, 0.0) - value).norm())
        .min(1.0 / norm)
}

#[cfg(test)]
mod tests {

    use crate::numerics::ApproxTolerance;

    #[test]
    fn new_preserves_caller_supplied_bounds() {
        let tolerance = ApproxTolerance::new(1.0e-8, 2.5e-7);

        assert_eq!(tolerance.absolute, 1.0e-8);
        assert_eq!(tolerance.relative, 2.5e-7);
    }

    #[test]
    fn educational_default_uses_documented_constants() {
        assert_eq!(
            ApproxTolerance::educational_default(),
            ApproxTolerance::new(1.0e-9, 1.0e-9)
        );
    }

    #[test]
    fn strict_is_tighter_than_the_educational_default() {
        let strict = ApproxTolerance::strict();
        let default = ApproxTolerance::educational_default();

        assert!(strict.absolute < default.absolute);
        assert!(strict.relative < default.relative);
    }

    #[test]
    fn loose_is_more_permissive_than_the_educational_default() {
        let loose = ApproxTolerance::loose();
        let default = ApproxTolerance::educational_default();

        assert!(loose.absolute > default.absolute);
        assert!(loose.relative > default.relative);
    }

    #[test]
    fn real_close_uses_absolute_tolerance_near_zero() {
        let tolerance = ApproxTolerance::new(1.0e-9, 1.0e-12);

        assert!(tolerance.real_close(0.0, 5.0e-10));
        assert!(!tolerance.real_close(0.0, 2.0e-9));
    }

    #[test]
    fn real_close_uses_relative_tolerance_at_large_scale() {
        let tolerance = ApproxTolerance::new(1.0e-12, 1.0e-6);

        assert!(tolerance.real_close(1_000_000.0, 1_000_000.5));
        assert!(!tolerance.real_close(1_000_000.0, 1_000_002.0));
    }
}
