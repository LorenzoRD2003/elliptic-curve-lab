use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    periods::{
        WeierstrassCubicRoots,
        legendre::{
            LegendreOrbitElementKind, LegendreParameter, LegendreParameterConditioning,
            LegendreReduction, conditioning::LegendreSingularityDiagnostics,
        },
    },
};
use crate::numerics::ApproxTolerance;

/// Structured report for one Legendre reduction.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreReductionReport {
    reduction: LegendreReduction,
    selected_orbit_element_relative_to_input_order: LegendreOrbitElementKind,
    conditioning: LegendreParameterConditioning,
    tolerance: ApproxTolerance,
    singularity_distance: f64,
}

impl LegendreReductionReport {
    /// Builds a Legendre reduction report from an already computed reduction.
    pub fn new(reduction: LegendreReduction, tolerance: ApproxTolerance) -> Self {
        let selected_orbit_element =
            reduction.selected_orbit_element_relative_to_input_order(tolerance);
        let diagnostics = LegendreSingularityDiagnostics::analyze(reduction.parameter(), tolerance);

        Self {
            reduction,
            selected_orbit_element_relative_to_input_order: selected_orbit_element,
            conditioning: diagnostics.conditioning(),
            tolerance,
            singularity_distance: diagnostics.singularity_distance(),
        }
    }

    /// Builds a report directly from a root triple.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        Ok(Self::new(
            LegendreReduction::from_roots(roots, tolerance)?,
            tolerance,
        ))
    }

    pub fn reduction(&self) -> &LegendreReduction {
        &self.reduction
    }

    pub fn parameter(&self) -> &LegendreParameter {
        self.reduction.parameter()
    }

    pub fn selected_orbit_element_relative_to_input_order(&self) -> LegendreOrbitElementKind {
        self.selected_orbit_element_relative_to_input_order
    }

    pub fn conditioning(&self) -> LegendreParameterConditioning {
        self.conditioning
    }

    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    pub fn singularity_distance(&self) -> f64 {
        self.singularity_distance
    }

    pub fn is_near_singular(&self) -> bool {
        self.conditioning.is_near_singular()
    }
}
