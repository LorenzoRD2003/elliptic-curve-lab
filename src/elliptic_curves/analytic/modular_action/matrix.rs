use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};

/// A validated matrix in the modular group `SL_2(ℤ)`.
///
/// This type stores an integer matrix
/// `γ = [[a, b], [c, d]]` with determinant `ad - bc = 1`.
/// It acts on the upper half-plane by fractional linear transformations
/// `γ(τ) = (aτ + b) / (cτ + d)`.
///
/// In particular, the standard generators satisfy
/// `S(τ) = -1/τ` and `T(τ) = τ + 1`.
///
/// The key geometric point is that this action does **not** usually keep the
/// complex number `τ` fixed. Instead, it changes the lattice basis used to
/// describe the same torus. If `Λ_τ = ℤ + ℤτ`, then
/// `Λ_{γ(τ)} = (1 / (cτ + d)) Λ_τ`, so the two lattices differ only by a
/// nonzero complex scaling factor. They therefore define isomorphic complex
/// tori. Equivalently, `τ` and `γ(τ)` represent the same point of the modular
/// quotient `SL_2(ℤ)\backslash\mathfrak H`.
///
/// The entries are machine integers, so this is an educational bounded model
/// of `SL_2(ℤ)` rather than an unbounded exact one.
///
/// Group operations use checked integer arithmetic. That means a mathematically
/// valid product or inverse can still fail with
/// [`AnalyticCurveError::InvalidModularMatrix`] if the intermediate `i128`
/// arithmetic would overflow.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModularMatrix {
    a: i128,
    b: i128,
    c: i128,
    d: i128,
}

impl ModularMatrix {
    /// Builds a validated modular matrix.
    ///
    /// The constructor accepts exactly the matrices in `SL_2(ℤ)`,
    /// meaning those with determinant `ad - bc = 1`.
    pub fn new(a: i128, b: i128, c: i128, d: i128) -> Result<Self, AnalyticCurveError> {
        let determinant = checked_determinant(a, b, c, d)?;
        if determinant != 1 {
            return Err(AnalyticCurveError::InvalidModularMatrix);
        }

        Ok(Self { a, b, c, d })
    }

    /// Returns the identity matrix.
    pub fn identity() -> Self {
        Self {
            a: 1,
            b: 0,
            c: 0,
            d: 1,
        }
    }

    /// Returns `S = [[0, -1], [1, 0]]`.
    ///
    /// Under the modular action, this generator satisfies `S(τ) = -1/τ`.
    pub fn s() -> Self {
        Self {
            a: 0,
            b: -1,
            c: 1,
            d: 0,
        }
    }

    /// Returns `T = [[1, 1], [0, 1]]`.
    ///
    /// Under the modular action, this generator satisfies `T(τ) = τ + 1`.
    pub fn t() -> Self {
        Self {
            a: 1,
            b: 1,
            c: 0,
            d: 1,
        }
    }

    /// Returns the top-left entry `a`.
    pub fn a(&self) -> i128 {
        self.a
    }

    /// Returns the top-right entry `b`.
    pub fn b(&self) -> i128 {
        self.b
    }

    /// Returns the bottom-left entry `c`.
    pub fn c(&self) -> i128 {
        self.c
    }

    /// Returns the bottom-right entry `d`.
    pub fn d(&self) -> i128 {
        self.d
    }

    /// Returns the determinant `ad - bc`.
    pub fn determinant(&self) -> i128 {
        checked_determinant(self.a, self.b, self.c, self.d)
            .expect("validated modular matrices keep determinant arithmetic in range")
    }

    /// Returns the product `self * other`.
    ///
    /// Matrix multiplication is performed with checked `i128` arithmetic.
    pub fn compose(&self, other: &Self) -> Result<Self, AnalyticCurveError> {
        let a = checked_add(checked_mul(self.a, other.a)?, checked_mul(self.b, other.c)?)?;
        let b = checked_add(checked_mul(self.a, other.b)?, checked_mul(self.b, other.d)?)?;
        let c = checked_add(checked_mul(self.c, other.a)?, checked_mul(self.d, other.c)?)?;
        let d = checked_add(checked_mul(self.c, other.b)?, checked_mul(self.d, other.d)?)?;

        Self::new(a, b, c, d)
    }

    /// Returns the inverse matrix `γ^{-1}`.
    ///
    /// For `γ = [[a, b], [c, d]]` in `SL_2(ℤ)`, the inverse is
    /// `[[d, -b], [-c, a]]`. The sign changes use checked negation.
    pub fn inverse(&self) -> Result<Self, AnalyticCurveError> {
        Self::new(self.d, checked_neg(self.b)?, checked_neg(self.c)?, self.a)
    }

    /// Applies the fractional linear transformation
    /// `τ ↦ (aτ + b) / (cτ + d)` to a point in the upper half-plane.
    ///
    /// For the standard generators this specializes to
    /// `S(τ) = -1/τ` and `T(τ) = τ + 1`.
    ///
    /// Intuitively, this produces a new coordinate for the same geometric
    /// complex torus, not a genuinely unrelated torus. So the output `γτ`
    /// should be read as “the same object written in a different modular
    /// coordinate”.
    ///
    /// For a valid matrix in `SL_2(ℤ)`, this action should preserve
    /// the upper half-plane mathematically. The implementation still validates
    /// the result numerically and reports a dedicated error if floating-point
    /// roundoff produces a non-positive imaginary part.
    pub fn apply(
        &self,
        tau: &UpperHalfPlanePoint,
    ) -> Result<UpperHalfPlanePoint, AnalyticCurveError> {
        let tau_value = *tau.tau();
        let numerator =
            Complex64::new(self.a as f64, 0.0) * tau_value + Complex64::new(self.b as f64, 0.0);
        let denominator =
            Complex64::new(self.c as f64, 0.0) * tau_value + Complex64::new(self.d as f64, 0.0);
        let denominator_norm_sqr = denominator.norm_sqr();

        if !denominator_norm_sqr.is_finite() || denominator_norm_sqr <= 0.0 {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        let image = numerator / denominator;
        if !image.re.is_finite() || !image.im.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        UpperHalfPlanePoint::new(image)
            .map_err(|_| AnalyticCurveError::NonPositiveImaginaryPartAfterModularAction)
    }
}

fn checked_mul(lhs: i128, rhs: i128) -> Result<i128, AnalyticCurveError> {
    lhs.checked_mul(rhs)
        .ok_or(AnalyticCurveError::InvalidModularMatrix)
}

fn checked_add(lhs: i128, rhs: i128) -> Result<i128, AnalyticCurveError> {
    lhs.checked_add(rhs)
        .ok_or(AnalyticCurveError::InvalidModularMatrix)
}

fn checked_sub(lhs: i128, rhs: i128) -> Result<i128, AnalyticCurveError> {
    lhs.checked_sub(rhs)
        .ok_or(AnalyticCurveError::InvalidModularMatrix)
}

fn checked_neg(value: i128) -> Result<i128, AnalyticCurveError> {
    value
        .checked_neg()
        .ok_or(AnalyticCurveError::InvalidModularMatrix)
}

fn checked_determinant(a: i128, b: i128, c: i128, d: i128) -> Result<i128, AnalyticCurveError> {
    checked_sub(checked_mul(a, d)?, checked_mul(b, c)?)
}
