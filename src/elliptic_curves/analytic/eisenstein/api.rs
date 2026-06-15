use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, LatticeSumTruncation,
    eisenstein::types::{EisensteinSumApprox, TruncationConvergenceReport},
};
use crate::numerics::ComplexDifferenceReport;

impl ComplexLattice {
    /// Approximates the lattice Eisenstein sum `G_k(Λ)` by a finite punctured
    /// square-box sum.
    pub fn eisenstein_sum(
        &self,
        weight: usize,
        truncation: LatticeSumTruncation,
    ) -> Result<EisensteinSumApprox, AnalyticCurveError> {
        if weight < 3 {
            return Err(AnalyticCurveError::InvalidEisensteinWeight);
        }

        let terms_used = truncation.nonzero_terms_in_square_box();

        if weight % 2 == 1 {
            return Ok(EisensteinSumApprox::new(
                weight,
                Complex64::new(0.0, 0.0),
                truncation,
                terms_used,
            ));
        }

        let value = self
            .nonzero_lattice_points_in_box(truncation.radius())
            .into_iter()
            .fold(Complex64::new(0.0, 0.0), |acc, point| {
                acc + Complex64::new(1.0, 0.0) / point.value.powu(weight as u32)
            });

        Ok(EisensteinSumApprox::new(
            weight, value, truncation, terms_used,
        ))
    }

    /// Approximates the classical weight-4 Eisenstein sum `G₄(Λ)`.
    pub fn g4_sum(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<EisensteinSumApprox, AnalyticCurveError> {
        self.eisenstein_sum(4, truncation)
    }

    /// Approximates the classical weight-6 Eisenstein sum `G₆(Λ)`.
    pub fn g6_sum(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<EisensteinSumApprox, AnalyticCurveError> {
        self.eisenstein_sum(6, truncation)
    }

    /// Compares two truncations of the same lattice Eisenstein sum.
    pub fn compare_eisenstein_truncations(
        &self,
        weight: usize,
        small: LatticeSumTruncation,
        large: LatticeSumTruncation,
    ) -> Result<TruncationConvergenceReport, AnalyticCurveError> {
        if large.radius() <= small.radius() {
            return Err(AnalyticCurveError::InvalidTruncationComparison);
        }

        let small_sum = self.eisenstein_sum(weight, small)?;
        let large_sum = self.eisenstein_sum(weight, large)?;
        Ok(TruncationConvergenceReport::new(
            small,
            large,
            ComplexDifferenceReport::new(*large_sum.value(), *small_sum.value()),
        ))
    }
}
