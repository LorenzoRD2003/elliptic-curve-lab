use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ApproxTolerance, ModularMatrix, UpperHalfPlanePoint,
    fundamental_domain::types::{
        FundamentalDomainReductionReport, FundamentalDomainReductionStatus,
        FundamentalDomainReductionStep, FundamentalDomainReductionStepReason,
    },
};
use crate::fields::complex_approx::ComplexApprox;
use num_bigint::BigInt;
use num_traits::FromPrimitive;

impl UpperHalfPlanePoint {
    /// Returns whether `self` lies in the classical standard fundamental
    /// domain `|Re(τ)| ≤ 1/2`, `|τ| ≥ 1` under the supplied numerical
    /// tolerance.
    pub fn is_in_standard_fundamental_domain(&self, tolerance: ApproxTolerance) -> bool {
        self.real_part().abs() <= 0.5 + tolerance.absolute
            && self.norm_sqr() >= 1.0 - tolerance.absolute
    }

    /// Reduces `self` toward the standard fundamental domain of `SL_2(ℤ)`.
    ///
    /// The routine alternates the classical translation and inversion moves,
    /// while recording the accumulated modular matrix and each concrete step.
    pub fn reduce_to_standard_fundamental_domain(
        &self,
        max_steps: usize,
    ) -> Result<FundamentalDomainReductionReport, AnalyticCurveError> {
        let original_tau = self.clone();
        let mut current = self.clone();
        let mut accumulated_matrix = ModularMatrix::identity();
        let mut steps = Vec::new();
        let tolerance = ComplexApprox::default_tolerance();

        loop {
            if current.is_in_standard_fundamental_domain(tolerance) {
                let status = if steps.is_empty() {
                    FundamentalDomainReductionStatus::AlreadyReduced
                } else {
                    FundamentalDomainReductionStatus::Reduced
                };

                return Ok(FundamentalDomainReductionReport::new(
                    original_tau,
                    current,
                    accumulated_matrix,
                    steps,
                    status,
                ));
            }

            if steps.len() >= max_steps {
                return Ok(FundamentalDomainReductionReport::new(
                    original_tau,
                    current,
                    accumulated_matrix,
                    steps,
                    FundamentalDomainReductionStatus::StepLimitReached,
                ));
            }

            let (step_matrix, reason) = next_reduction_step(&current, tolerance)?;
            let before = current.clone();
            current = step_matrix.apply(&current)?;
            accumulated_matrix = step_matrix.compose(&accumulated_matrix)?;
            steps.push(FundamentalDomainReductionStep::new(
                step_matrix,
                before,
                current.clone(),
                reason,
            ));
        }
    }
}

fn next_reduction_step(
    tau: &UpperHalfPlanePoint,
    tolerance: ApproxTolerance,
) -> Result<(ModularMatrix, FundamentalDomainReductionStepReason), AnalyticCurveError> {
    if let Some(shift) = centered_strip_shift(tau, tolerance) {
        return Ok((
            ModularMatrix::new(1, -shift, 0, 1)?,
            FundamentalDomainReductionStepReason::RealPartOutsideCenteredStrip,
        ));
    }

    if tau_norm_is_strictly_less_than_one(tau, tolerance) {
        return Ok((
            ModularMatrix::s(),
            FundamentalDomainReductionStepReason::NormLessThanOne,
        ));
    }

    Err(AnalyticCurveError::NumericalComparisonFailed)
}

fn centered_strip_shift(tau: &UpperHalfPlanePoint, tolerance: ApproxTolerance) -> Option<BigInt> {
    let real_part = tau.real_part();
    if real_part.abs() <= 0.5 + tolerance.absolute {
        return None;
    }

    BigInt::from_f64(real_part.round())
}

fn tau_norm_is_strictly_less_than_one(
    tau: &UpperHalfPlanePoint,
    tolerance: ApproxTolerance,
) -> bool {
    tau.norm_sqr() < 1.0 - tolerance.absolute
}
