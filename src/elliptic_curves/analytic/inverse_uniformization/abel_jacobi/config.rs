use crate::ApproxTolerance;

/// Contour-family policy for the current Abel-Jacobi integral.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LegendreContourStrategy {
    /// Use the current deterministic Legendre-side `segment + ray` contour.
    CanonicalSegmentThenRay,
}

/// Validation-truncation policy for the final torus-to-curve roundtrip check.
///
/// This is intentionally separate from the actual Abel-Jacobi quadrature
/// budget, so callers can experiment with forward-validation sensitivity
/// without implicitly changing the inverse integral itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AbelJacobiValidationPolicy {
    pub lattice_truncation_radius: usize,
    pub function_truncation_radius: usize,
}

impl AbelJacobiValidationPolicy {
    /// Returns the baseline validation policy for educational experiments.
    pub fn educational_default() -> Self {
        Self {
            lattice_truncation_radius: 16,
            function_truncation_radius: 14,
        }
    }

    /// Returns a tighter validation policy for more delicate runs.
    pub fn strict() -> Self {
        Self {
            lattice_truncation_radius: 24,
            function_truncation_radius: 22,
        }
    }

    /// Returns a lighter validation policy for coarse exploratory work.
    pub fn loose() -> Self {
        Self {
            lattice_truncation_radius: 12,
            function_truncation_radius: 10,
        }
    }
}

/// Configuration for the pedagogical Abel-Jacobi inverse map `(x, y) -> z in C / Λ`.
///
/// Current scope:
/// - one explicit contour-strategy selector for the improper integral from `∞` to `x`
/// - numerical branch-continuation corrections for the square root
/// - one independent validation policy for the final torus-to-curve roundtrip
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbelJacobiConfig {
    /// Tolerance policy for branch tracking and final validation.
    pub tolerance: ApproxTolerance,
    /// Contour-family policy used by the Legendre-side integral stage.
    pub legendre_contour_strategy: LegendreContourStrategy,
    /// Total Simpson-quadrature budget for the `segment + ray` integral.
    pub integration_steps: usize,
    /// Number of samples used on the finite segment when scoring candidate
    /// contours against the branch locus.
    pub segment_samples: usize,
    /// Number of samples used on the compactified ray when scoring candidate
    /// contours against the branch locus.
    pub ray_samples: usize,
    /// Maximum number of sign-flip corrections allowed while continuing the
    /// square-root branch along the contour.
    pub max_branch_adjustments: usize,
    /// Maximum number of lattice-side correction attempts allowed while
    /// reducing the recovered complex value modulo the period lattice.
    pub max_lattice_corrections: usize,
    /// Truncation policy used by the final roundtrip validation through
    /// `(wp, wp')`.
    pub validation_policy: AbelJacobiValidationPolicy,
}

impl AbelJacobiConfig {
    /// Returns the baseline preset for educational experiments.
    pub fn educational_default() -> Self {
        Self {
            tolerance: ApproxTolerance::educational_default(),
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            integration_steps: 256,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 8,
            max_lattice_corrections: 4,
            validation_policy: AbelJacobiValidationPolicy::educational_default(),
        }
    }

    /// Returns a tighter preset for more delicate inverse-uniformization runs.
    pub fn strict() -> Self {
        Self {
            tolerance: ApproxTolerance::strict(),
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            integration_steps: 512,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 16,
            max_lattice_corrections: 8,
            validation_policy: AbelJacobiValidationPolicy::strict(),
        }
    }

    /// Returns a more permissive preset for coarse exploratory work.
    pub fn loose() -> Self {
        Self {
            tolerance: ApproxTolerance::loose(),
            legendre_contour_strategy: LegendreContourStrategy::CanonicalSegmentThenRay,
            integration_steps: 128,
            segment_samples: 32,
            ray_samples: 32,
            max_branch_adjustments: 4,
            max_lattice_corrections: 2,
            validation_policy: AbelJacobiValidationPolicy::loose(),
        }
    }
}

/// Terminal status for one Abel-Jacobi point-recovery attempt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbelJacobiRecoveryStatus {
    Succeeded,
    HitBranchAdjustmentLimit,
    HitLatticeCorrectionLimit,
    ValidationFailed,
    Failed,
}

/// Numerical metadata for one Abel-Jacobi inverse-uniformization run.
#[derive(Clone, Debug, PartialEq)]
pub struct AbelJacobiRecoveryMetadata {
    status: AbelJacobiRecoveryStatus,
    integration_steps_used: usize,
    branch_adjustments_used: usize,
    lattice_corrections_used: usize,
    tolerance: ApproxTolerance,
    validation_x_residual_norm: Option<f64>,
    validation_y_residual_norm: Option<f64>,
}

impl AbelJacobiRecoveryMetadata {
    /// Builds one explicit Abel-Jacobi numerical metadata bundle.
    pub fn new(
        status: AbelJacobiRecoveryStatus,
        integration_steps_used: usize,
        branch_adjustments_used: usize,
        lattice_corrections_used: usize,
        tolerance: ApproxTolerance,
        validation_x_residual_norm: Option<f64>,
        validation_y_residual_norm: Option<f64>,
    ) -> Self {
        Self {
            status,
            integration_steps_used,
            branch_adjustments_used,
            lattice_corrections_used,
            tolerance,
            validation_x_residual_norm,
            validation_y_residual_norm,
        }
    }

    /// Returns the interpreted terminal status.
    pub fn status(&self) -> AbelJacobiRecoveryStatus {
        self.status
    }

    /// Returns the quadrature step count used by the current run.
    pub fn integration_steps_used(&self) -> usize {
        self.integration_steps_used
    }

    /// Returns how many branch-continuation corrections were applied.
    pub fn branch_adjustments_used(&self) -> usize {
        self.branch_adjustments_used
    }

    /// Returns how many lattice-side corrections were applied before the
    /// result was accepted.
    pub fn lattice_corrections_used(&self) -> usize {
        self.lattice_corrections_used
    }

    /// Returns the tolerance policy used by the current run.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the validation residual for the recovered `x = wp(z)` value
    /// when such a validation step was carried out.
    pub fn validation_x_residual_norm(&self) -> Option<f64> {
        self.validation_x_residual_norm
    }

    /// Returns the validation residual for the recovered `y = wp'(z)` value
    /// when such a validation step was carried out.
    pub fn validation_y_residual_norm(&self) -> Option<f64> {
        self.validation_y_residual_norm
    }

    /// Returns whether the run ended in the success state.
    pub fn succeeded(&self) -> bool {
        self.status == AbelJacobiRecoveryStatus::Succeeded
    }
}
