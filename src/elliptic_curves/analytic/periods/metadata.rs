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
    cardano_product_residual_norm: Option<f64>,
    cardano_discriminant: Option<Complex64>,
    selected_u_branch_index: Option<usize>,
    selected_v_branch_index: Option<usize>,
}

impl NumericalRecoveryMetadata {
    /// Builds one numerical-recovery metadata bundle from explicit fields.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
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
            cardano_product_residual_norm: None,
            cardano_discriminant: None,
            selected_u_branch_index: None,
            selected_v_branch_index: None,
        }
    }

    /// Attaches Cardano-branch diagnostics gathered during cubic-root
    /// recovery.
    pub fn with_cardano_diagnostics(
        mut self,
        cardano_discriminant: Complex64,
        cardano_product_residual_norm: f64,
        selected_u_branch_index: usize,
        selected_v_branch_index: usize,
    ) -> Self {
        self.cardano_discriminant = Some(cardano_discriminant);
        self.cardano_product_residual_norm = Some(cardano_product_residual_norm);
        self.selected_u_branch_index = Some(selected_u_branch_index);
        self.selected_v_branch_index = Some(selected_v_branch_index);
        self
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

    /// Returns the residual norm of the selected Cardano branch condition
    /// `uv ≈ -p/3`, when the run recorded it.
    pub fn cardano_product_residual_norm(&self) -> Option<f64> {
        self.cardano_product_residual_norm
    }

    /// Returns the Cardano discriminant
    /// `(q/2)^2 + (p/3)^3`, when the run recorded it.
    pub fn cardano_discriminant(&self) -> Option<&Complex64> {
        self.cardano_discriminant.as_ref()
    }

    /// Returns the selected branch index for `u`, when the run recorded it.
    pub fn selected_u_branch_index(&self) -> Option<usize> {
        self.selected_u_branch_index
    }

    /// Returns the selected branch index for `v`, when the run recorded it.
    pub fn selected_v_branch_index(&self) -> Option<usize> {
        self.selected_v_branch_index
    }

    /// Returns whether the selected Cardano pair used the principal
    /// cube-root branch for both `u` and `v`, when the run recorded the
    /// branch indices.
    pub fn used_principal_cardano_branches(&self) -> Option<bool> {
        match (self.selected_u_branch_index, self.selected_v_branch_index) {
            (Some(u_index), Some(v_index)) => Some(u_index == 0 && v_index == 0),
            _ => None,
        }
    }

    /// Returns whether the run ended in the succeeded state.
    pub fn succeeded(&self) -> bool {
        self.status == PeriodRecoveryStatus::Succeeded
    }
}
