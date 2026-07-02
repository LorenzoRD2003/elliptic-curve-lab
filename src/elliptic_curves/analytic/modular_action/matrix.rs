use num_bigint::BigInt;
use num_complex::Complex64;
use num_traits::{One, ToPrimitive, Zero};

use crate::elliptic_curves::analytic::{AnalyticCurveError, UpperHalfPlanePoint};

/// A validated matrix in the modular group `SL₂(ℤ)`.
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
/// quotient `SL₂(ℤ)\backslash\mathfrak H`.
///
/// Entries are exact unbounded integers. The only lossy step is applying the
/// matrix to an approximate complex number, where entries must be representable
/// as finite `f64` values.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModularMatrix {
    a: BigInt,
    b: BigInt,
    c: BigInt,
    d: BigInt,
}

impl ModularMatrix {
    /// Builds a validated modular matrix.
    ///
    /// The constructor accepts exactly the matrices in `SL₂(ℤ)`, meaning those
    /// with determinant `ad - bc = 1`.
    pub fn new(
        a: impl Into<BigInt>,
        b: impl Into<BigInt>,
        c: impl Into<BigInt>,
        d: impl Into<BigInt>,
    ) -> Result<Self, AnalyticCurveError> {
        let matrix = Self {
            a: a.into(),
            b: b.into(),
            c: c.into(),
            d: d.into(),
        };

        if matrix.determinant() != BigInt::one() {
            return Err(AnalyticCurveError::InvalidModularMatrix);
        }

        Ok(matrix)
    }

    /// Returns the identity matrix.
    pub fn identity() -> Self {
        Self {
            a: BigInt::one(),
            b: BigInt::zero(),
            c: BigInt::zero(),
            d: BigInt::one(),
        }
    }

    /// Returns `S = [[0, -1], [1, 0]]`.
    ///
    /// Under the modular action, this generator satisfies `S(τ) = -1/τ`.
    pub fn s() -> Self {
        Self {
            a: BigInt::zero(),
            b: BigInt::from(-1),
            c: BigInt::one(),
            d: BigInt::zero(),
        }
    }

    /// Returns `T = [[1, 1], [0, 1]]`.
    ///
    /// Under the modular action, this generator satisfies `T(τ) = τ + 1`.
    pub fn t() -> Self {
        Self {
            a: BigInt::one(),
            b: BigInt::one(),
            c: BigInt::zero(),
            d: BigInt::one(),
        }
    }

    /// Returns the top-left entry `a`.
    pub fn a(&self) -> &BigInt {
        &self.a
    }

    /// Returns the top-right entry `b`.
    pub fn b(&self) -> &BigInt {
        &self.b
    }

    /// Returns the bottom-left entry `c`.
    pub fn c(&self) -> &BigInt {
        &self.c
    }

    /// Returns the bottom-right entry `d`.
    pub fn d(&self) -> &BigInt {
        &self.d
    }

    /// Returns the determinant `ad - bc`.
    pub fn determinant(&self) -> BigInt {
        &self.a * &self.d - &self.b * &self.c
    }

    /// Returns the product `self * other`.
    pub fn compose(&self, other: &Self) -> Result<Self, AnalyticCurveError> {
        Self::new(
            &self.a * &other.a + &self.b * &other.c,
            &self.a * &other.b + &self.b * &other.d,
            &self.c * &other.a + &self.d * &other.c,
            &self.c * &other.b + &self.d * &other.d,
        )
    }

    /// Returns the inverse matrix `γ⁻¹`.
    ///
    /// For `γ = [[a, b], [c, d]]` in `SL₂(ℤ)`, the inverse is
    /// `[[d, -b], [-c, a]]`.
    pub fn inverse(&self) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.d.clone(),
            -self.b.clone(),
            -self.c.clone(),
            self.a.clone(),
        )
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
    /// For a valid matrix in `SL₂(ℤ)`, this action should preserve the upper
    /// half-plane mathematically. The implementation still validates the
    /// result numerically and reports a dedicated error if floating-point
    /// conversion or roundoff leaves the represented chart.
    pub fn apply(
        &self,
        tau: &UpperHalfPlanePoint,
    ) -> Result<UpperHalfPlanePoint, AnalyticCurveError> {
        let a = big_int_to_f64(&self.a)?;
        let b = big_int_to_f64(&self.b)?;
        let c = big_int_to_f64(&self.c)?;
        let d = big_int_to_f64(&self.d)?;

        let tau_value = *tau.tau();
        let numerator = Complex64::new(a, 0.0) * tau_value + Complex64::new(b, 0.0);
        let denominator = Complex64::new(c, 0.0) * tau_value + Complex64::new(d, 0.0);
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

fn big_int_to_f64(value: &BigInt) -> Result<f64, AnalyticCurveError> {
    value
        .to_f64()
        .filter(|value| value.is_finite())
        .ok_or(AnalyticCurveError::NumericalComparisonFailed)
}
