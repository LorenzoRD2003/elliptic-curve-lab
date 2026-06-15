use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ComplexLattice, EllipticFunctionTruncation,
    LatticeSumTruncation, TorusToCurveMapResult, UpperHalfPlanePoint,
    lattice::HasAnalyticLatticeContext,
};
use crate::numerics::ApproxTolerance;

/// Aggregated experiment for the analytic uniformization map
/// `ℂ / Λ → E(ℂ)`, `z ↦ (℘(z), ℘′(z))`.
///
/// This report keeps one lattice/curve pair together with several sampled
/// torus representatives and the corresponding mapped curve points.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct UniformizationExperimentReport {
    tau: UpperHalfPlanePoint,
    lattice: ComplexLattice,
    curve: AnalyticWeierstrassCurve,
    sampled_points: Vec<TorusToCurveMapResult>,
    all_points_lie_on_curve: bool,
}

#[allow(dead_code)]
impl UniformizationExperimentReport {
    /// Builds one uniformization experiment from explicit complex
    /// representatives in `ℂ`.
    ///
    /// The sampled points may include genuine lattice points; those are mapped
    /// to [`super::super::AnalyticCurvePoint::Infinity`] through the same pole
    /// convention already used by [`ComplexLattice::map_torus_point_to_curve`].
    pub(crate) fn from_sample_points(
        tau: UpperHalfPlanePoint,
        sample_points: Vec<Complex64>,
        invariant_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        let lattice = ComplexLattice::from_tau(tau.clone());
        let curve = AnalyticWeierstrassCurve::from_lattice(&lattice, invariant_truncation)?;
        let sampled_points = sample_points
            .into_iter()
            .map(|z| {
                lattice.map_torus_point_to_curve(
                    z,
                    invariant_truncation,
                    function_truncation,
                    tolerance,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let all_points_lie_on_curve = sampled_points.iter().all(|point| point.lies_on_curve());

        Ok(Self {
            tau,
            lattice,
            curve,
            sampled_points,
            all_points_lie_on_curve,
        })
    }

    /// Returns the original upper-half-plane parameter.
    pub(crate) fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the underlying lattice.
    pub(crate) fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the common analytic cubic used in the experiment.
    pub(crate) fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the sampled torus-to-curve map results.
    pub(crate) fn sampled_points(&self) -> &[TorusToCurveMapResult] {
        &self.sampled_points
    }

    /// Returns whether every sampled point was accepted as lying on the curve.
    pub(crate) fn all_points_lie_on_curve(&self) -> bool {
        self.all_points_lie_on_curve
    }
}

impl HasAnalyticLatticeContext for UniformizationExperimentReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}
