use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, ComplexLattice, EllipticFunctionTruncation,
    LatticeSumTruncation, torsion::AnalyticTorsionPointApprox,
};

impl ComplexLattice {
    /// Maps every torus `n`-torsion point to the analytic Weierstrass curve
    /// attached to the same lattice.
    ///
    /// The torus side is enumerated in the same lexicographic `(a, b)` order as
    /// [`Self::torus_n_torsion_points`]. Each representative
    /// `z = (a/n)ω₁ + (b/n)ω₂ (mod Λ)` is then sent through the analytic
    /// correspondence `z ↦ (℘(z), ℘′(z))`.
    ///
    /// The identity class `(a, b) = (0, 0)` is treated specially: because it
    /// represents a lattice point, `℘` and `℘′` have poles there, so the map lands
    /// at [`AnalyticCurvePoint::Infinity`] instead of returning an evaluation error.
    ///
    /// Complexity: `Θ(n²(r_inv² + r_fun²))`, where `r_inv` is the invariant
    /// truncation radius and `r_fun` is the elliptic-function truncation radius.
    pub fn map_torus_torsion_to_curve(
        &self,
        n: usize,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<Vec<AnalyticTorsionPointApprox>, AnalyticCurveError> {
        self.torus_n_torsion_points(n)?
            .into_iter()
            .map(|torus_point| {
                let map_result = self.map_torus_point_to_curve(
                    *torus_point.z(),
                    invariant_truncation,
                    function_truncation,
                    tolerance,
                )?;

                Ok(AnalyticTorsionPointApprox::new(
                    torus_point,
                    map_result.point().clone(),
                    map_result.membership_report().clone(),
                ))
            })
            .collect()
    }

    /// Maps the primitive torus `n`-torsion points to the analytic Weierstrass
    /// curve attached to the same lattice.
    ///
    /// This uses the same torus-side ordering and the same analytic correspondence
    /// as [`Self::map_torus_torsion_to_curve`], but keeps only the primitive classes,
    /// meaning those with exact torus order `n` or equivalently
    /// `gcd(a, b, n) = 1`.
    ///
    /// The identity class is included only in the special case `n = 1`, where it
    /// is primitive by the current convention and still maps to
    /// [`AnalyticCurvePoint::Infinity`].
    ///
    /// Complexity: `Θ(n²(r_inv² + r_fun²))`, because the current implementation
    /// maps the full torus `n`-torsion grid first and then filters it.
    pub fn map_primitive_torus_torsion_to_curve(
        &self,
        n: usize,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<Vec<AnalyticTorsionPointApprox>, AnalyticCurveError> {
        Ok(self
            .map_torus_torsion_to_curve(n, invariant_truncation, function_truncation, tolerance)?
            .into_iter()
            .filter(|point| point.index().is_primitive())
            .collect())
    }
}
