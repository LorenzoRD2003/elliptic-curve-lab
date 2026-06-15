use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError,
    periods::{
        WeierstrassCubicRoots,
        legendre::{
            LegendreParameterConditioning, LegendreParameterOrbit,
            conditioning::LegendreSingularityDiagnostics,
        },
    },
};
use crate::numerics::ApproxTolerance;

/// One chosen Legendre parameter `λ` for a cubic normalized to
/// `y² = x(x - 1)(x - λ)`.
#[derive(Clone, Debug, PartialEq)]
pub struct LegendreParameter {
    lambda: Complex64,
}

impl LegendreParameter {
    /// Builds a finite Legendre parameter from an already chosen `λ`.
    pub fn new(lambda: Complex64) -> Result<Self, AnalyticCurveError> {
        if !lambda.is_finite()
            || lambda == Complex64::new(0.0, 0.0)
            || lambda == Complex64::new(1.0, 0.0)
        {
            return Err(AnalyticCurveError::InvalidLegendreModulus);
        }

        Ok(Self { lambda })
    }

    /// Builds one deterministic Legendre parameter from an unordered cubic
    /// root triple.
    pub fn from_roots(
        roots: &WeierstrassCubicRoots,
        tolerance: ApproxTolerance,
    ) -> Result<Self, AnalyticCurveError> {
        roots.validate_distinct(tolerance)?;
        let selected = roots.choose_legendre_candidate(tolerance)?;
        selected.into_parameter(tolerance)
    }

    pub fn lambda(&self) -> &Complex64 {
        &self.lambda
    }

    pub fn orbit(&self) -> LegendreParameterOrbit {
        LegendreParameterOrbit::from_parameter(self)
    }

    pub fn one_minus_lambda(&self) -> Complex64 {
        Complex64::new(1.0, 0.0) - self.lambda
    }

    pub fn is_near_zero(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_zero()
    }

    pub fn is_near_one(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_one()
    }

    pub fn is_near_singular(&self, tolerance: ApproxTolerance) -> bool {
        LegendreSingularityDiagnostics::analyze(self, tolerance).is_near_singular()
    }

    pub fn conditioning(&self, tolerance: ApproxTolerance) -> LegendreParameterConditioning {
        LegendreSingularityDiagnostics::analyze(self, tolerance).conditioning()
    }
}
