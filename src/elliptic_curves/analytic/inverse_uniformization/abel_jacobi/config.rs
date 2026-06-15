use crate::elliptic_curves::analytic::{
    AnalyticCurveError, EllipticFunctionTruncation, LatticeSumTruncation, PeriodRecoveryConfig,
};
use crate::numerics::ApproxTolerance;

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
///
/// The stored radii are private and validated. Each radius must be positive,
/// because the corresponding lattice or elliptic-function truncation would
/// otherwise collapse to a degenerate educational regime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AbelJacobiValidationPolicy {
    lattice_truncation_radius: usize,
    function_truncation_radius: usize,
}

impl AbelJacobiValidationPolicy {
    /// Builds an explicit validation policy from positive truncation radii.
    pub fn new(
        lattice_truncation_radius: usize,
        function_truncation_radius: usize,
    ) -> Result<Self, AnalyticCurveError> {
        LatticeSumTruncation::new(lattice_truncation_radius)?;
        EllipticFunctionTruncation::new(function_truncation_radius)?;

        Ok(Self {
            lattice_truncation_radius,
            function_truncation_radius,
        })
    }

    /// Returns the baseline validation policy for educational experiments.
    pub fn educational_default() -> Self {
        Self::new(16, 14).expect("educational Abel-Jacobi validation policy is valid")
    }

    /// Returns a tighter validation policy for more delicate runs.
    pub fn strict() -> Self {
        Self::new(24, 22).expect("strict Abel-Jacobi validation policy is valid")
    }

    /// Returns a lighter validation policy for coarse exploratory work.
    pub fn loose() -> Self {
        Self::new(12, 10).expect("loose Abel-Jacobi validation policy is valid")
    }

    /// Returns the lattice-sum truncation radius used by the final roundtrip
    /// validation.
    pub fn lattice_truncation_radius(&self) -> usize {
        self.lattice_truncation_radius
    }

    /// Returns the elliptic-function truncation radius used by the final
    /// roundtrip validation.
    pub fn function_truncation_radius(&self) -> usize {
        self.function_truncation_radius
    }

    /// Returns a copy of this policy with a different positive lattice
    /// truncation radius.
    pub fn with_lattice_truncation_radius(
        self,
        lattice_truncation_radius: usize,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(lattice_truncation_radius, self.function_truncation_radius)
    }

    /// Returns a copy of this policy with a different positive
    /// elliptic-function truncation radius.
    pub fn with_function_truncation_radius(
        self,
        function_truncation_radius: usize,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(self.lattice_truncation_radius, function_truncation_radius)
    }

    /// Returns the lattice-sum truncation used by the final roundtrip
    /// validation.
    pub(crate) fn lattice_truncation(self) -> Result<LatticeSumTruncation, AnalyticCurveError> {
        LatticeSumTruncation::new(self.lattice_truncation_radius)
    }

    /// Returns the elliptic-function truncation used by the final roundtrip
    /// validation.
    pub(crate) fn function_truncation(
        self,
    ) -> Result<EllipticFunctionTruncation, AnalyticCurveError> {
        EllipticFunctionTruncation::new(self.function_truncation_radius)
    }
}

/// Configuration for the pedagogical Abel-Jacobi inverse map `(x, y) -> z in C / Λ`.
///
/// Current scope:
/// - one explicit contour-strategy selector for the improper integral from `∞` to `x`
/// - numerical branch-continuation corrections for the square root
/// - one independent validation policy for the final torus-to-curve roundtrip
///
/// All stored numerical budgets are validated to be positive. In normal use,
/// callers should start from `educational_default()`, `strict()`, or `loose()`
/// and then refine individual knobs with the `with_*` helpers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbelJacobiConfig {
    tolerance: ApproxTolerance,
    legendre_contour_strategy: LegendreContourStrategy,
    integration_steps: usize,
    segment_samples: usize,
    ray_samples: usize,
    max_branch_adjustments: usize,
    max_lattice_corrections: usize,
    validation_policy: AbelJacobiValidationPolicy,
}

impl AbelJacobiConfig {
    /// Builds an explicit Abel-Jacobi configuration from validated numerical
    /// knobs.
    pub fn new(
        tolerance: ApproxTolerance,
        legendre_contour_strategy: LegendreContourStrategy,
        integration_steps: usize,
        segment_samples: usize,
        ray_samples: usize,
        max_branch_adjustments: usize,
        max_lattice_corrections: usize,
        validation_policy: AbelJacobiValidationPolicy,
    ) -> Result<Self, AnalyticCurveError> {
        for value in [
            integration_steps,
            segment_samples,
            ray_samples,
            max_branch_adjustments,
            max_lattice_corrections,
        ] {
            if value == 0 {
                return Err(AnalyticCurveError::InvalidAbelJacobiConfig);
            }
        }

        Ok(Self {
            tolerance,
            legendre_contour_strategy,
            integration_steps,
            segment_samples,
            ray_samples,
            max_branch_adjustments,
            max_lattice_corrections,
            validation_policy,
        })
    }

    /// Returns the baseline preset for educational experiments.
    pub fn educational_default() -> Self {
        Self::new(
            ApproxTolerance::educational_default(),
            LegendreContourStrategy::CanonicalSegmentThenRay,
            256,
            32,
            32,
            8,
            4,
            AbelJacobiValidationPolicy::educational_default(),
        )
        .expect("educational Abel-Jacobi config is valid")
    }

    /// Returns a tighter preset for more delicate inverse-uniformization runs.
    pub fn strict() -> Self {
        Self::new(
            ApproxTolerance::strict(),
            LegendreContourStrategy::CanonicalSegmentThenRay,
            512,
            32,
            32,
            16,
            8,
            AbelJacobiValidationPolicy::strict(),
        )
        .expect("strict Abel-Jacobi config is valid")
    }

    /// Returns a more permissive preset for coarse exploratory work.
    pub fn loose() -> Self {
        Self::new(
            ApproxTolerance::loose(),
            LegendreContourStrategy::CanonicalSegmentThenRay,
            128,
            32,
            32,
            4,
            2,
            AbelJacobiValidationPolicy::loose(),
        )
        .expect("loose Abel-Jacobi config is valid")
    }

    /// Returns the tolerance policy used for branch tracking and final
    /// validation.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the contour-family policy used by the Legendre-side integral
    /// stage.
    pub fn legendre_contour_strategy(&self) -> LegendreContourStrategy {
        self.legendre_contour_strategy
    }

    /// Returns the total Simpson-quadrature budget for the `segment + ray`
    /// integral.
    pub fn integration_steps(&self) -> usize {
        self.integration_steps
    }

    /// Returns the number of finite-segment samples used while scoring
    /// candidate contours.
    pub fn segment_samples(&self) -> usize {
        self.segment_samples
    }

    /// Returns the number of compactified-ray samples used while scoring
    /// candidate contours.
    pub fn ray_samples(&self) -> usize {
        self.ray_samples
    }

    /// Returns the maximum number of sign-flip corrections allowed during
    /// square-root branch continuation.
    pub fn max_branch_adjustments(&self) -> usize {
        self.max_branch_adjustments
    }

    /// Returns the maximum number of lattice-side correction attempts allowed
    /// while reducing the recovered complex value modulo the period lattice.
    pub fn max_lattice_corrections(&self) -> usize {
        self.max_lattice_corrections
    }

    /// Returns the validation policy used by the final roundtrip check
    /// through `(wp, wp')`.
    pub fn validation_policy(&self) -> AbelJacobiValidationPolicy {
        self.validation_policy
    }

    /// Returns a copy of this config with a different tolerance policy.
    pub fn with_tolerance(self, tolerance: ApproxTolerance) -> Result<Self, AnalyticCurveError> {
        Self::new(
            tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different contour strategy.
    pub fn with_legendre_contour_strategy(
        self,
        legendre_contour_strategy: LegendreContourStrategy,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different positive quadrature
    /// budget.
    pub fn with_integration_steps(
        self,
        integration_steps: usize,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            integration_steps,
            self.segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different positive segment-sample
    /// count.
    pub fn with_segment_samples(self, segment_samples: usize) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different positive ray-sample
    /// count.
    pub fn with_ray_samples(self, ray_samples: usize) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different positive branch-adjustment
    /// budget.
    pub fn with_max_branch_adjustments(
        self,
        max_branch_adjustments: usize,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            self.ray_samples,
            max_branch_adjustments,
            self.max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different positive
    /// lattice-correction budget.
    pub fn with_max_lattice_corrections(
        self,
        max_lattice_corrections: usize,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            max_lattice_corrections,
            self.validation_policy,
        )
    }

    /// Returns a copy of this config with a different validation policy.
    pub fn with_validation_policy(
        self,
        validation_policy: AbelJacobiValidationPolicy,
    ) -> Result<Self, AnalyticCurveError> {
        Self::new(
            self.tolerance,
            self.legendre_contour_strategy,
            self.integration_steps,
            self.segment_samples,
            self.ray_samples,
            self.max_branch_adjustments,
            self.max_lattice_corrections,
            validation_policy,
        )
    }

    /// Derives the period-recovery configuration used internally by the
    /// current Abel-Jacobi implementation when it needs cubic-root recovery.
    pub(crate) fn derived_root_recovery_config(self) -> PeriodRecoveryConfig {
        PeriodRecoveryConfig::new(
            self.tolerance,
            (self.max_branch_adjustments * 2).max(8),
            6,
            self.integration_steps.max(8),
            self.max_lattice_corrections.max(1),
            8,
        )
        .expect("derived Abel-Jacobi root-recovery config must stay valid")
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
    pub(crate) fn new(
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
