use num_complex::Complex64;

use crate::fields::ComplexApprox;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, UpperHalfPlanePoint, lattice::ComplexLattice,
};

impl ComplexLattice {
    /// Validates that two complex basis vectors define a positively oriented,
    /// non-degenerate lattice basis.
    ///
    /// The current convention accepts exactly those pairs with strictly
    /// positive oriented area, equivalently those with `Im(ω₂ / ω₁) > 0`
    /// when `ω₁ != 0`.
    ///
    /// If the oriented area is numerically too close to zero, this returns
    /// [`AnalyticCurveError::DegenerateLattice`]. If the basis spans a genuine
    /// parallelogram but with non-positive orientation, this returns
    /// [`AnalyticCurveError::NonPositiveLatticeOrientation`].
    pub fn validate_basis(omega1: Complex64, omega2: Complex64) -> Result<(), AnalyticCurveError> {
        let area = Self::oriented_area_from_basis(omega1, omega2);
        let tolerance = ComplexApprox::default_tolerance().absolute;

        if area > tolerance {
            Ok(())
        } else if area < -tolerance {
            Err(AnalyticCurveError::NonPositiveLatticeOrientation)
        } else {
            Err(AnalyticCurveError::DegenerateLattice)
        }
    }

    /// Builds a validated complex lattice from an ordered basis.
    pub fn new(omega1: Complex64, omega2: Complex64) -> Result<Self, AnalyticCurveError> {
        Self::validate_basis(omega1, omega2)?;

        Ok(Self { omega1, omega2 })
    }

    /// Builds the standard lattice `ℤ + ℤτ` from a point `τ` in the upper
    /// half-plane.
    pub fn from_tau(tau: UpperHalfPlanePoint) -> Self {
        Self::new(Complex64::new(1.0, 0.0), *tau.tau())
            .expect("a valid upper-half-plane parameter defines a positive basis")
    }

    /// Returns the first lattice basis vector `ω₁`.
    pub fn omega1(&self) -> &Complex64 {
        &self.omega1
    }

    /// Returns the second lattice basis vector `ω₂`.
    pub fn omega2(&self) -> &Complex64 {
        &self.omega2
    }

    /// Returns the associated upper-half-plane parameter `τ = ω₂ / ω₁`.
    ///
    /// Because construction already enforces a positively oriented
    /// non-degenerate basis, this conversion should succeed for valid lattice
    /// values. It still returns a `Result` to keep the API honest and local to
    /// the analytic error surface.
    pub fn tau(&self) -> Result<UpperHalfPlanePoint, AnalyticCurveError> {
        UpperHalfPlanePoint::new(self.omega2 / self.omega1)
    }

    /// Returns the oriented real area of the fundamental parallelogram.
    ///
    /// Concretely, if `ω₁ = a + bi` and `ω₂ = c + di`, this is the determinant
    /// `ad - bc`. Positive values correspond to the accepted orientation
    /// convention used by [`Self::new`], zero corresponds to a degenerate
    /// basis, and negative values correspond to the rejected opposite
    /// orientation.
    pub fn oriented_area(&self) -> f64 {
        Self::oriented_area_from_basis(self.omega1, self.omega2)
    }

    /// Returns the covolume of the lattice.
    ///
    /// The covolume is the geometric area of the fundamental parallelogram,
    /// so it is the absolute value of the oriented area. Under the current
    /// positive-orientation convention this agrees numerically with
    /// [`Self::oriented_area`], but the distinction is still mathematically
    /// useful and worth documenting explicitly.
    pub fn covolume(&self) -> f64 {
        self.oriented_area().abs()
    }

    pub(super) fn oriented_area_from_basis(omega1: Complex64, omega2: Complex64) -> f64 {
        omega1.re * omega2.im - omega1.im * omega2.re
    }
}
