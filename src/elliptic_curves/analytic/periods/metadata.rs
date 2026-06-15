use num_complex::Complex64;

use crate::numerics::ApproxTolerance;

/// High-level strategy used by one numerical period-recovery attempt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodRecoveryMethod {
    /// Recover periods through a Legendre-normalized AGM route.
    AgmViaLegendre,
    /// Recover periods through direct numerical path integration.
    NumericalPathIntegral,
    /// Combine more than one numerical route in one recovery attempt.
    Hybrid,
}

/// Outcome status for one numerical period-recovery attempt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeriodRecoveryStatus {
    /// The requested recovery workflow completed successfully.
    Succeeded,
    /// The workflow stopped because one iteration budget was exhausted.
    HitIterationLimit,
    /// The workflow could not choose a consistent branch representative.
    BranchChoiceAmbiguous,
    /// The workflow produced candidate periods that failed validation checks.
    ValidationFailed,
    /// The workflow failed for another reason not represented above.
    Failed,
}

/// Structured execution metadata for one numerical period-recovery run.
///
/// This stores the resolved recovery route, its outcome, the work counters for
/// the main numerical phases, and one optional validation residual norm when a
/// final comparison was available.
#[derive(Clone, Debug, PartialEq)]
pub struct NumericalRecoveryMetadata {
    resolved_method: PeriodRecoveryMethod,
    status: PeriodRecoveryStatus,
    newton_iterations_used: usize,
    agm_iterations_used: usize,
    integration_steps_used: usize,
    branch_lattice_searches_used: usize,
    tolerance: ApproxTolerance,
    validation_residual_norm: Option<f64>,
}

impl NumericalRecoveryMetadata {
    /// Builds one numerical-recovery metadata bundle from explicit fields.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        resolved_method: PeriodRecoveryMethod,
        status: PeriodRecoveryStatus,
        newton_iterations_used: usize,
        agm_iterations_used: usize,
        integration_steps_used: usize,
        branch_lattice_searches_used: usize,
        tolerance: ApproxTolerance,
        validation_residual_norm: Option<f64>,
    ) -> Self {
        Self {
            resolved_method,
            status,
            newton_iterations_used,
            agm_iterations_used,
            integration_steps_used,
            branch_lattice_searches_used,
            tolerance,
            validation_residual_norm,
        }
    }

    /// Returns the resolved numerical recovery route used for this run.
    pub fn resolved_method(&self) -> PeriodRecoveryMethod {
        self.resolved_method
    }

    /// Returns the final status of this recovery run.
    pub fn status(&self) -> PeriodRecoveryStatus {
        self.status
    }

    /// Returns the number of Newton iterations consumed.
    pub fn newton_iterations_used(&self) -> usize {
        self.newton_iterations_used
    }

    /// Returns the number of AGM iterations consumed.
    pub fn agm_iterations_used(&self) -> usize {
        self.agm_iterations_used
    }

    /// Returns the number of numerical integration steps consumed.
    pub fn integration_steps_used(&self) -> usize {
        self.integration_steps_used
    }

    /// Returns the number of nearby lattice-branch searches consumed.
    pub fn branch_lattice_searches_used(&self) -> usize {
        self.branch_lattice_searches_used
    }

    /// Returns the tolerance policy used by the run.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the final validation residual norm, when one was recorded.
    pub fn validation_residual_norm(&self) -> Option<f64> {
        self.validation_residual_norm
    }

    /// Returns whether the run ended in the succeeded state.
    pub fn succeeded(&self) -> bool {
        self.status == PeriodRecoveryStatus::Succeeded
    }

    pub(crate) fn from_root_and_integral_reports(
        root_metadata: &NumericalRecoveryMetadata,
        integral_report: &crate::elliptic_curves::analytic::LegendrePeriodIntegralReport,
        recovered_tau: Complex64,
    ) -> Self {
        let agm_iterations_used = integral_report.k_lambda().metadata().agm_iterations_used()
            + integral_report
                .k_complementary()
                .metadata()
                .agm_iterations_used();
        let status = [
            root_metadata.status(),
            integral_report.k_lambda().metadata().status(),
            integral_report.k_complementary().metadata().status(),
        ]
        .into_iter()
        .find(|status| *status != PeriodRecoveryStatus::Succeeded)
        .unwrap_or(PeriodRecoveryStatus::Succeeded);
        let validation_residual_norm = Some(
            root_metadata
                .validation_residual_norm()
                .unwrap_or(0.0)
                .max((recovered_tau - *integral_report.tau_candidate()).norm()),
        );

        NumericalRecoveryMetadata::new(
            PeriodRecoveryMethod::Hybrid,
            status,
            root_metadata.newton_iterations_used(),
            agm_iterations_used,
            0,
            root_metadata.branch_lattice_searches_used(),
            root_metadata.tolerance(),
            validation_residual_norm,
        )
    }
}
