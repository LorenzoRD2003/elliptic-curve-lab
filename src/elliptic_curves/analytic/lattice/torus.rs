use num_complex::Complex64;

use super::{
    super::AnalyticCurveError, ComplexLattice, ComplexTorusPoint,
    FundamentalParallelogramCoordinate,
};

impl ComplexTorusPoint {
    pub(super) fn new(coordinate: FundamentalParallelogramCoordinate) -> Self {
        Self { coordinate }
    }

    /// Returns the canonical reduced coordinate inside `[0, 1) × [0, 1)`.
    pub fn coordinate(&self) -> &FundamentalParallelogramCoordinate {
        &self.coordinate
    }
}

impl ComplexLattice {
    /// Compares two torus points relative to this ambient lattice.
    ///
    /// Equality on `ℂ / Λ` is meaningful only after choosing the quotient
    /// lattice `Λ`. For that reason [`ComplexTorusPoint`] deliberately does
    /// not implement `PartialEq`. This method provides the explicit
    /// lattice-relative comparison surface instead.
    ///
    /// Under the current canonical-coordinate representation, two torus
    /// points are equal exactly when their reduced fundamental-parallelogram
    /// coordinates agree.
    pub fn torus_points_eq(&self, lhs: &ComplexTorusPoint, rhs: &ComplexTorusPoint) -> bool {
        let _ = self;
        lhs.coordinate() == rhs.coordinate()
    }

    /// Converts fundamental coordinates into the corresponding complex point.
    ///
    /// For a lattice basis `ω₁, ω₂`, the coordinate pair `(u, v)` is sent to
    /// `uω₁ + vω₂`. Because `FundamentalParallelogramCoordinate` validates its
    /// inputs, the returned complex number is always the chosen
    /// representative inside the half-open fundamental parallelogram.
    pub fn point_from_fundamental_coordinates(
        &self,
        coord: FundamentalParallelogramCoordinate,
    ) -> Complex64 {
        self.omega1 * coord.u() + self.omega2 * coord.v()
    }

    /// Reduces a complex point to canonical coordinates in `[0, 1) × [0, 1)`.
    ///
    /// Every `z ∈ ℂ` has a unique class in the quotient `ℂ / Λ`. Relative to
    /// the ordered basis `ω₁, ω₂`, we first solve for the real coefficients
    /// `u, v` in
    ///
    /// `z = uω₁ + vω₂`,
    ///
    /// then reduce `u` and `v` modulo `1` to land in the chosen half-open
    /// fundamental parallelogram.
    ///
    /// This method returns an error only when floating-point arithmetic
    /// produces non-finite intermediate coordinates or the final normalized
    /// pair cannot be confirmed to lie in the canonical region.
    pub fn reduce_complex_point_to_fundamental_coordinates(
        &self,
        z: Complex64,
    ) -> Result<FundamentalParallelogramCoordinate, AnalyticCurveError> {
        if !z.re.is_finite() || !z.im.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        let determinant = self.oriented_area();
        let raw_u = (z.re * self.omega2.im - z.im * self.omega2.re) / determinant;
        let raw_v = (self.omega1.re * z.im - self.omega1.im * z.re) / determinant;

        if !raw_u.is_finite() || !raw_v.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        FundamentalParallelogramCoordinate::reduce_mod_unit_square(raw_u, raw_v)
    }

    /// Reduces a complex point to its canonical torus representative.
    ///
    /// This is the quotient-space view of
    /// [`Self::reduce_complex_point_to_fundamental_coordinates`]: instead of
    /// returning just the reduced pair `(u, v)`, it packages that pair as a
    /// canonical point of the torus `ℂ / Λ`.
    pub fn reduce_complex_point_to_torus_point(
        &self,
        z: Complex64,
    ) -> Result<ComplexTorusPoint, AnalyticCurveError> {
        self.reduce_complex_point_to_fundamental_coordinates(z)
            .map(ComplexTorusPoint::new)
    }
}
