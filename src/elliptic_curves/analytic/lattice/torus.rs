use num_complex::Complex64;

use super::{
    super::{AnalyticCurveError, ApproxTolerance, ComplexApproxComparison},
    ComplexLattice, ComplexModuloLatticeComparison, ComplexTorusPoint,
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

    /// Compares two complex representatives approximately modulo this lattice.
    ///
    /// The current implementation searches over the square box
    /// `-search_radius ≤ m, n ≤ search_radius` and chooses the lattice shift
    /// `mω₁ + nω₂` minimizing
    ///
    /// `|z_left - (z_right + mω₁ + nω₂)|`.
    ///
    /// The returned report stores that best shift together with the resulting
    /// [`ComplexApproxComparison`].
    ///
    /// Complexity: `Θ(search_radius²)`, because the searched lattice box is
    /// traversed explicitly.
    pub fn compare_complex_points_mod_lattice_approx(
        &self,
        z_left: Complex64,
        z_right: Complex64,
        search_radius: usize,
        tolerance: ApproxTolerance,
    ) -> Result<ComplexModuloLatticeComparison, AnalyticCurveError> {
        if !z_left.is_finite() || !z_right.is_finite() {
            return Err(AnalyticCurveError::NumericalComparisonFailed);
        }

        let best_shift = self
            .lattice_points_in_box(search_radius)
            .into_iter()
            .min_by(|lhs, rhs| {
                let lhs_norm = (z_left - (z_right + lhs.value)).norm();
                let rhs_norm = (z_left - (z_right + rhs.value)).norm();
                lhs_norm.total_cmp(&rhs_norm)
            })
            .expect("lattice_points_in_box always includes the origin");

        Ok(ComplexModuloLatticeComparison {
            original_right: z_right,
            best_shift: best_shift.clone(),
            comparison: ComplexApproxComparison::new(z_left, z_right + best_shift.value, tolerance),
            search_radius,
        })
    }
}

impl ComplexModuloLatticeComparison {
    /// Returns the left-hand representative being compared.
    pub fn left_representative(&self) -> &Complex64 {
        self.comparison.left()
    }

    /// Returns the original unshifted right-hand representative.
    pub fn right_representative(&self) -> &Complex64 {
        &self.original_right
    }

    /// Returns the best lattice shift found within the searched box.
    pub fn best_shift(&self) -> &super::LatticeIndexPoint {
        &self.best_shift
    }

    /// Returns the right-hand representative after applying the best shift.
    pub fn shifted_right_representative(&self) -> &Complex64 {
        self.comparison.right()
    }

    /// Returns the final comparison payload after the best shift was applied.
    pub fn comparison(&self) -> &ComplexApproxComparison {
        &self.comparison
    }

    /// Returns the searched square-box radius in lattice coordinates.
    pub fn search_radius(&self) -> usize {
        self.search_radius
    }

    /// Returns the residual before any lattice shift was applied.
    pub fn unshifted_difference(&self) -> Complex64 {
        *self.left_representative() - self.original_right
    }

    /// Returns the Euclidean norm of the unshifted residual.
    pub fn unshifted_difference_norm(&self) -> f64 {
        self.unshifted_difference().norm()
    }

    /// Returns the shifted residual `z_left - (z_right + λ_best)`.
    pub fn shifted_difference(&self) -> &Complex64 {
        self.comparison.difference()
    }

    /// Returns the Euclidean norm of the shifted residual.
    pub fn shifted_difference_norm(&self) -> f64 {
        self.comparison.absolute_difference()
    }

    /// Returns whether the best shifted comparison passed the supplied
    /// tolerance.
    pub fn agrees_approximately(&self) -> bool {
        self.comparison.agrees_approximately()
    }
}
