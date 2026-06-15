use num_complex::Complex64;

use crate::elliptic_curves::analytic::periods::LegendreParameter;
use crate::fields::complex_approx::ComplexApprox;
use crate::numerics::{
    ApproxTolerance, projective_unit_singularity_distance, reciprocal_singularity_threshold,
};

/// Coarse numerical conditioning class for one chosen Legendre parameter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendreParameterConditioning {
    /// The chosen `λ` is not close to `0`, `1`, or `∞` under the tolerance.
    Generic,
    /// The chosen `λ` is numerically close to `0`.
    NearZero,
    /// The chosen `λ` is numerically close to `1`.
    NearOne,
    /// The chosen `λ` is numerically close to `∞`, detected through `1/λ`.
    NearInfinity,
}

impl LegendreParameterConditioning {
    /// Returns whether this conditioning class lies near the singular
    /// Legendre locus `{0, 1, ∞}`.
    pub fn is_near_singular(self) -> bool {
        self != Self::Generic
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LegendreSingularityDiagnostics {
    conditioning: LegendreParameterConditioning,
    singularity_distance: f64,
    is_near_zero: bool,
    is_near_one: bool,
}

impl LegendreSingularityDiagnostics {
    pub(crate) fn analyze(parameter: &LegendreParameter, tolerance: ApproxTolerance) -> Self {
        let is_near_zero = ComplexApprox::eq_with_tolerance(
            parameter.lambda(),
            &Complex64::new(0.0, 0.0),
            tolerance,
        );
        let is_near_one = ComplexApprox::eq_with_tolerance(
            parameter.lambda(),
            &Complex64::new(1.0, 0.0),
            tolerance,
        );
        let is_near_infinity =
            parameter.lambda().norm() >= reciprocal_singularity_threshold(tolerance);
        let conditioning = if is_near_zero {
            LegendreParameterConditioning::NearZero
        } else if is_near_one {
            LegendreParameterConditioning::NearOne
        } else if is_near_infinity {
            LegendreParameterConditioning::NearInfinity
        } else {
            LegendreParameterConditioning::Generic
        };

        Self {
            conditioning,
            singularity_distance: projective_unit_singularity_distance(parameter.lambda()),
            is_near_zero,
            is_near_one,
        }
    }

    pub(crate) fn conditioning(self) -> LegendreParameterConditioning {
        self.conditioning
    }

    pub(crate) fn singularity_distance(self) -> f64 {
        self.singularity_distance
    }

    pub(crate) fn is_near_zero(self) -> bool {
        self.is_near_zero
    }

    pub(crate) fn is_near_one(self) -> bool {
        self.is_near_one
    }

    pub(crate) fn is_near_singular(self) -> bool {
        self.conditioning.is_near_singular()
    }
}
