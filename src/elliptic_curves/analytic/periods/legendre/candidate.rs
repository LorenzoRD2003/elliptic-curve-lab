use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, periods::LegendreParameter};
use crate::numerics::{ApproxTolerance, projective_unit_singularity_distance};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LegendreCandidate {
    permutation: [usize; 3],
    lambda: Complex64,
    singular_distance: f64,
}

impl LegendreCandidate {
    pub(crate) fn new(permutation: [usize; 3], lambda: Complex64) -> Self {
        Self {
            permutation,
            lambda,
            singular_distance: projective_unit_singularity_distance(&lambda),
        }
    }

    pub(crate) fn permutation(self) -> [usize; 3] {
        self.permutation
    }

    pub(crate) fn is_better_than(
        self,
        current_best: Option<Self>,
        tolerance: ApproxTolerance,
    ) -> bool {
        let Some(best) = current_best else {
            return true;
        };

        if !tolerance.real_close(self.singular_distance, best.singular_distance) {
            return self.singular_distance > best.singular_distance;
        }

        let candidate_norm = self.lambda.norm();
        let best_norm = best.lambda.norm();
        if !tolerance.real_close(candidate_norm, best_norm) {
            return candidate_norm < best_norm;
        }

        if !tolerance.real_close(self.lambda.re, best.lambda.re) {
            return self.lambda.re < best.lambda.re;
        }

        self.lambda.im < best.lambda.im
    }

    pub(crate) fn into_parameter(
        self,
        tolerance: ApproxTolerance,
    ) -> Result<LegendreParameter, AnalyticCurveError> {
        let parameter = LegendreParameter::new(self.lambda)?;
        if parameter.is_near_singular(tolerance) {
            return Err(AnalyticCurveError::InvalidLegendreModulus);
        }

        Ok(parameter)
    }
}
