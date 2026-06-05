use super::{AnalyticCurveError, ApproxTolerance, ModularMatrix, UpperHalfPlanePoint};
use crate::fields::ComplexApprox;

/// Records one actual modular transformation used while reducing `τ` toward
/// the standard fundamental domain.
#[derive(Clone, Debug, PartialEq)]
pub struct FundamentalDomainReductionStep {
    applied_matrix: ModularMatrix,
    before: UpperHalfPlanePoint,
    after: UpperHalfPlanePoint,
    reason: FundamentalDomainReductionStepReason,
}

impl FundamentalDomainReductionStep {
    /// Returns the matrix applied at this step.
    pub fn applied_matrix(&self) -> ModularMatrix {
        self.applied_matrix
    }

    /// Returns the upper-half-plane point before the step.
    pub fn before(&self) -> &UpperHalfPlanePoint {
        &self.before
    }

    /// Returns the upper-half-plane point after the step.
    pub fn after(&self) -> &UpperHalfPlanePoint {
        &self.after
    }

    /// Returns why this modular transformation was applied.
    pub fn reason(&self) -> FundamentalDomainReductionStepReason {
        self.reason
    }
}

/// Explains why one modular step was applied during reduction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FundamentalDomainReductionStepReason {
    /// A translation `τ ↦ τ - n` was needed because `Re(τ)` lay outside the
    /// centered strip `|Re(τ)| ≤ 1/2`.
    RealPartOutsideCenteredStrip,
    /// An inversion `τ ↦ -1/τ` was needed because `|τ| < 1`.
    NormLessThanOne,
}

/// Final status of one attempted reduction to the standard fundamental domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FundamentalDomainReductionStatus {
    /// The input already satisfied the standard domain conditions, so no
    /// modular step was needed.
    AlreadyReduced,
    /// The algorithm reached a point satisfying the current domain test.
    Reduced,
    /// The algorithm stopped after `max_steps` transformations without yet
    /// reaching the target region.
    StepLimitReached,
}

/// Structured report for reducing `τ` to the classical standard fundamental
/// domain of `SL_2(ℤ)`.
///
/// The target region is `|Re(τ)| ≤ 1/2`, `|τ| ≥ 1` with a numerical boundary
/// convention: values within the default [`ComplexApprox`] tolerance of the
/// boundary are treated as already inside.
///
/// Geometrically, this reduction chooses a more canonical representative for
/// the modular class of `τ` in the quotient `SL_2(ℤ) \ ℍ`.
#[derive(Clone, Debug, PartialEq)]
pub struct FundamentalDomainReductionReport {
    original_tau: UpperHalfPlanePoint,
    reduced_tau: UpperHalfPlanePoint,
    accumulated_matrix: ModularMatrix,
    steps: Vec<FundamentalDomainReductionStep>,
    status: FundamentalDomainReductionStatus,
}

impl FundamentalDomainReductionReport {
    /// Returns the original upper-half-plane input.
    pub fn original_tau(&self) -> &UpperHalfPlanePoint {
        &self.original_tau
    }

    /// Returns the final representative reached by the reduction routine.
    pub fn reduced_tau(&self) -> &UpperHalfPlanePoint {
        &self.reduced_tau
    }

    /// Returns the accumulated modular matrix `γ` such that the final point is
    /// `γ(τ)` for the original input `τ`.
    pub fn accumulated_matrix(&self) -> ModularMatrix {
        self.accumulated_matrix
    }

    /// Returns the sequence of actual modular transformations that were
    /// applied.
    pub fn steps(&self) -> &[FundamentalDomainReductionStep] {
        &self.steps
    }

    /// Returns the terminal status of the reduction attempt.
    pub fn status(&self) -> FundamentalDomainReductionStatus {
        self.status
    }

    /// Returns whether the final point was accepted as lying in the target
    /// standard fundamental domain.
    pub fn is_reduced(&self) -> bool {
        matches!(
            self.status,
            FundamentalDomainReductionStatus::AlreadyReduced
                | FundamentalDomainReductionStatus::Reduced
        )
    }
}

/// Reduces one point `τ ∈ ℍ` toward the standard fundamental domain of
/// `SL_2(ℤ)`.
///
/// The routine alternates two classical moves:
///
/// - translate by an integer until `Re(τ)` lies in the centered strip
///   `|Re(τ)| ≤ 1/2`
/// - if necessary, apply `S(τ) = -1/τ` to move from `|τ| < 1` to `|τ| > 1`
///
/// The result keeps a complete step history together with the accumulated
/// matrix mapping the original input to the final representative.
pub fn reduce_tau_to_standard_fundamental_domain(
    tau: UpperHalfPlanePoint,
    max_steps: usize,
) -> Result<FundamentalDomainReductionReport, AnalyticCurveError> {
    let mut current = tau.clone();
    let mut accumulated_matrix = ModularMatrix::identity();
    let mut steps = Vec::new();
    let tolerance = ComplexApprox::default_tolerance();

    loop {
        if is_in_standard_fundamental_domain(&current, tolerance) {
            let status = if steps.is_empty() {
                FundamentalDomainReductionStatus::AlreadyReduced
            } else {
                FundamentalDomainReductionStatus::Reduced
            };

            return Ok(FundamentalDomainReductionReport {
                original_tau: tau,
                reduced_tau: current,
                accumulated_matrix,
                steps,
                status,
            });
        }

        if steps.len() >= max_steps {
            return Ok(FundamentalDomainReductionReport {
                original_tau: tau,
                reduced_tau: current,
                accumulated_matrix,
                steps,
                status: FundamentalDomainReductionStatus::StepLimitReached,
            });
        }

        let (step_matrix, reason) = next_reduction_step(&current, tolerance)?;
        let before = current.clone();
        current = step_matrix.apply(&current)?;
        accumulated_matrix = step_matrix.compose(&accumulated_matrix)?;
        steps.push(FundamentalDomainReductionStep {
            applied_matrix: step_matrix,
            before,
            after: current.clone(),
            reason,
        });
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

/// Returns whether `τ` lies in the classical standard fundamental domain
/// `|Re(τ)| ≤ 1/2`, `|τ| ≥ 1` under the supplied numerical tolerance.
pub fn is_in_standard_fundamental_domain(
    tau: &UpperHalfPlanePoint,
    tolerance: ApproxTolerance,
) -> bool {
    tau.real_part().abs() <= 0.5 + tolerance.absolute && tau.norm_sqr() >= 1.0 - tolerance.absolute
}

fn centered_strip_shift(tau: &UpperHalfPlanePoint, tolerance: ApproxTolerance) -> Option<i128> {
    let real_part = tau.real_part();
    if real_part.abs() <= 0.5 + tolerance.absolute {
        return None;
    }

    Some(real_part.round() as i128)
}

fn tau_norm_is_strictly_less_than_one(
    tau: &UpperHalfPlanePoint,
    tolerance: ApproxTolerance,
) -> bool {
    tau.norm_sqr() < 1.0 - tolerance.absolute
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::{
        FundamentalDomainReductionStatus, FundamentalDomainReductionStepReason,
        is_in_standard_fundamental_domain, reduce_tau_to_standard_fundamental_domain,
    };
    use crate::elliptic_curves::analytic::{ApproxTolerance, ModularMatrix, UpperHalfPlanePoint};

    #[test]
    fn already_reduced_tau_reports_no_steps() {
        let tau = UpperHalfPlanePoint::tau_i();
        let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), 8)
            .expect("reduction should work");

        assert_eq!(report.original_tau(), &tau);
        assert_eq!(report.reduced_tau(), &tau);
        assert_eq!(report.accumulated_matrix(), ModularMatrix::identity());
        assert!(report.steps().is_empty());
        assert_eq!(
            report.status(),
            FundamentalDomainReductionStatus::AlreadyReduced
        );
        assert!(report.is_reduced());
    }

    #[test]
    fn translation_step_moves_tau_into_centered_strip() {
        let tau = UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap();
        let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), 8)
            .expect("reduction should work");

        assert_eq!(report.original_tau(), &tau);
        assert_eq!(report.status(), FundamentalDomainReductionStatus::Reduced);
        assert!(!report.steps().is_empty());
        assert_eq!(
            report.steps()[0].reason(),
            FundamentalDomainReductionStepReason::RealPartOutsideCenteredStrip
        );
        assert!(report.reduced_tau().real_part().abs() <= 0.5 + 1.0e-12);
        assert!(report.is_reduced());
    }

    #[test]
    fn inversion_step_applies_s_when_norm_is_less_than_one() {
        let tau = UpperHalfPlanePoint::from_re_im(0.1, 0.3).unwrap();
        let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), 8)
            .expect("reduction should work");

        assert!(!report.steps().is_empty());
        assert_eq!(
            report.steps()[0].reason(),
            FundamentalDomainReductionStepReason::NormLessThanOne
        );
        assert!(report.reduced_tau().norm_sqr() >= 1.0 - 1.0e-12);
        assert!(report.is_reduced());
    }

    #[test]
    fn step_limit_reached_is_reported_honestly() {
        let tau = UpperHalfPlanePoint::from_re_im(1.2, 1.0).unwrap();
        let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), 0)
            .expect("reduction should work");

        assert_eq!(report.original_tau(), &tau);
        assert_eq!(report.reduced_tau(), &tau);
        assert_eq!(report.accumulated_matrix(), ModularMatrix::identity());
        assert!(report.steps().is_empty());
        assert_eq!(
            report.status(),
            FundamentalDomainReductionStatus::StepLimitReached
        );
        assert!(!report.is_reduced());
    }

    #[test]
    fn accumulated_matrix_maps_original_tau_to_reduced_tau() {
        let tau = UpperHalfPlanePoint::from_re_im(1.2, 0.3).unwrap();
        let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), 8)
            .expect("reduction should work");

        let image = report
            .accumulated_matrix()
            .apply(&tau)
            .expect("accumulated matrix should act on tau");

        assert!((image.real_part() - report.reduced_tau().real_part()).abs() < 1.0e-12);
        assert!((image.imaginary_part() - report.reduced_tau().imaginary_part()).abs() < 1.0e-12);
    }

    #[test]
    fn explicit_domain_membership_helper_accepts_classical_points() {
        assert!(is_in_standard_fundamental_domain(
            &UpperHalfPlanePoint::tau_i(),
            ApproxTolerance::strict(),
        ));
        assert!(is_in_standard_fundamental_domain(
            &UpperHalfPlanePoint::tau_rho(),
            ApproxTolerance::strict(),
        ));
    }

    #[test]
    fn explicit_domain_membership_helper_respects_tolerance_near_the_boundary() {
        let tau = UpperHalfPlanePoint::from_re_im(0.5005, 1.0).unwrap();

        assert!(!is_in_standard_fundamental_domain(
            &tau,
            ApproxTolerance::strict(),
        ));
        assert!(is_in_standard_fundamental_domain(
            &tau,
            ApproxTolerance::new(1.0e-3, 1.0e-3),
        ));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(40))]

        #[test]
        fn generic_upper_half_plane_points_either_reduce_or_hit_the_requested_limit(
            re in -4.0f64..4.0,
            im in 0.1f64..4.0,
            max_steps in 0usize..8,
        ) {
            let tau = UpperHalfPlanePoint::from_re_im(re, im).unwrap();
            let report = reduce_tau_to_standard_fundamental_domain(tau.clone(), max_steps).unwrap();

            prop_assert_eq!(report.original_tau(), &tau);
            prop_assert!(report.steps().len() <= max_steps);

            match report.status() {
                FundamentalDomainReductionStatus::AlreadyReduced => {
                    prop_assert!(report.steps().is_empty());
                    prop_assert_eq!(report.reduced_tau(), &tau);
                }
                FundamentalDomainReductionStatus::Reduced => {
                    prop_assert!(report.reduced_tau().real_part().abs() <= 0.5 + 1.0e-9);
                    prop_assert!(report.reduced_tau().norm_sqr() >= 1.0 - 1.0e-9);
                }
                FundamentalDomainReductionStatus::StepLimitReached => {
                    prop_assert_eq!(report.steps().len(), max_steps);
                }
            }
        }
    }
}
