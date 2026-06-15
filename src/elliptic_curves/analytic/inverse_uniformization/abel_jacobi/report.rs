use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AbelJacobiRecoveryMetadata, AnalyticCurveError, AnalyticCurvePoint, AnalyticWeierstrassCurve,
    ComplexTorusPoint, PeriodBasisRecoveryReport, RecoveredPeriodBasis, UpperHalfPlanePoint,
    inverse_uniformization::abel_jacobi::{
        AbelJacobiRoundtripValidationReport, LegendreContourStrategy,
    },
};
use crate::numerics::{ApproxTolerance, ComplexApproxComparison};

/// Initial square-root branch choice used to start the Abel-Jacobi contour
/// integration in Legendre coordinates.
///
/// The integrand contains `sqrt(X(X-1)(X-λ))`, so the numerical routine must
/// choose one sign before branch continuation by continuity can begin.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbelJacobiInitialBranchChoice {
    Principal,
    Alternate,
}

/// Decomposes the numerical Abel-Jacobi contour integral into the pieces that
/// are actually accumulated by the current algorithm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbelJacobiIntegralDecomposition {
    initial_branch_choice: AbelJacobiInitialBranchChoice,
    segment_integral: Complex64,
    ray_integral: Complex64,
    tail_correction: Complex64,
}

impl AbelJacobiIntegralDecomposition {
    /// Builds one explicit Abel-Jacobi integral decomposition bundle.
    pub(crate) fn new(
        initial_branch_choice: AbelJacobiInitialBranchChoice,
        segment_integral: Complex64,
        ray_integral: Complex64,
        tail_correction: Complex64,
    ) -> Result<Self, AnalyticCurveError> {
        if !segment_integral.is_finite()
            || !ray_integral.is_finite()
            || !tail_correction.is_finite()
        {
            return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
        }

        Ok(Self {
            initial_branch_choice,
            segment_integral,
            ray_integral,
            tail_correction,
        })
    }

    /// Returns the initial square-root branch selected before continuity
    /// tracking along the contour begins.
    pub fn initial_branch_choice(&self) -> AbelJacobiInitialBranchChoice {
        self.initial_branch_choice
    }

    /// Returns the contribution of the finite segment from `X` to the anchor.
    pub fn segment_integral(&self) -> &Complex64 {
        &self.segment_integral
    }

    /// Returns the contribution of the compactified outgoing ray before the
    /// asymptotic tail correction.
    pub fn ray_integral(&self) -> &Complex64 {
        &self.ray_integral
    }

    /// Returns the asymptotic tail correction added beyond the sampled ray.
    pub fn tail_correction(&self) -> &Complex64 {
        &self.tail_correction
    }
}

/// Groups the numerical bookkeeping attached to one Abel-Jacobi integral
/// approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbelJacobiIntegralNumerics {
    integration_steps_used: usize,
    branch_adjustments_used: usize,
    tolerance: ApproxTolerance,
}

impl AbelJacobiIntegralNumerics {
    /// Builds one explicit Abel-Jacobi numerical-bookkeeping bundle.
    pub(crate) fn new(
        integration_steps_used: usize,
        branch_adjustments_used: usize,
        tolerance: ApproxTolerance,
    ) -> Self {
        Self {
            integration_steps_used,
            branch_adjustments_used,
            tolerance,
        }
    }

    /// Returns the quadrature step count used by the current run.
    pub fn integration_steps_used(&self) -> usize {
        self.integration_steps_used
    }

    /// Returns how many branch-continuation corrections were applied.
    pub fn branch_adjustments_used(&self) -> usize {
        self.branch_adjustments_used
    }

    /// Returns the comparison tolerance used by the current run.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }
}

/// Deterministic contour data used by the current Abel-Jacobi quadrature.
///
/// The present implementation uses a `segment + ray` contour in the Legendre
/// `X`-plane. This report records the concrete geometric choices so callers
/// can inspect which representative path was actually used.
#[derive(Clone, Debug, PartialEq)]
pub struct AbelJacobiContourReport {
    legendre_contour_strategy: LegendreContourStrategy,
    start: Complex64,
    anchor: Complex64,
    theta: f64,
    radius: f64,
    tail_length: f64,
    min_distance_to_branch_points: f64,
}

impl AbelJacobiContourReport {
    /// Builds one explicit Abel-Jacobi contour report.
    pub(crate) fn new(
        legendre_contour_strategy: LegendreContourStrategy,
        start: Complex64,
        anchor: Complex64,
        theta: f64,
        radius: f64,
        tail_length: f64,
        min_distance_to_branch_points: f64,
    ) -> Result<Self, AnalyticCurveError> {
        if !start.is_finite()
            || !anchor.is_finite()
            || !theta.is_finite()
            || !radius.is_finite()
            || !tail_length.is_finite()
            || !min_distance_to_branch_points.is_finite()
            || radius < 0.0
            || tail_length < 0.0
            || min_distance_to_branch_points < 0.0
        {
            return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
        }

        Ok(Self {
            legendre_contour_strategy,
            start,
            anchor,
            theta,
            radius,
            tail_length,
            min_distance_to_branch_points,
        })
    }

    /// Returns the contour-family strategy used to choose this report's
    /// concrete path.
    pub fn legendre_contour_strategy(&self) -> LegendreContourStrategy {
        self.legendre_contour_strategy
    }

    /// Returns the starting Legendre coordinate `X`.
    pub fn start(&self) -> &Complex64 {
        &self.start
    }

    /// Returns the contour anchor point.
    pub fn anchor(&self) -> &Complex64 {
        &self.anchor
    }

    /// Returns the chosen ray angle in radians.
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Returns the anchor radius `R`.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the radial length sampled on the outgoing ray before the
    /// asymptotic tail correction is applied.
    pub fn tail_length(&self) -> f64 {
        self.tail_length
    }

    /// Returns the minimum sampled distance from the contour to the branch
    /// locus `{0, 1, λ}`.
    pub fn min_distance_to_branch_points(&self) -> f64 {
        self.min_distance_to_branch_points
    }
}

/// Approximate value of the Abel-Jacobi integral before quotient reduction
/// modulo the recovered period lattice.
///
/// This stage intentionally stops short of interpreting the result as a torus
/// point. It records only the approximate complex integral
/// `z = ∫_x^∞ dt / √(4 t^3 - g₂ t - g₃)` together with the numerical
/// work spent on the quadrature and branch-continuation layers.
#[derive(Clone, Debug, PartialEq)]
pub struct AbelJacobiIntegralApprox {
    curve: AnalyticWeierstrassCurve,
    point: AnalyticCurvePoint,
    contour: AbelJacobiContourReport,
    value: Complex64,
    decomposition: AbelJacobiIntegralDecomposition,
    numerics: AbelJacobiIntegralNumerics,
}

impl AbelJacobiIntegralApprox {
    /// Builds one explicit Abel-Jacobi integral approximation report.
    pub(crate) fn new(
        curve: AnalyticWeierstrassCurve,
        point: AnalyticCurvePoint,
        contour: AbelJacobiContourReport,
        value: Complex64,
        decomposition: AbelJacobiIntegralDecomposition,
        numerics: AbelJacobiIntegralNumerics,
    ) -> Result<Self, AnalyticCurveError> {
        if !value.is_finite() {
            return Err(AnalyticCurveError::AbelJacobiIntegrationFailed);
        }

        Ok(Self {
            curve,
            point,
            contour,
            value,
            decomposition,
            numerics,
        })
    }

    /// Returns the source analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the source curve point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        &self.point
    }

    /// Returns the contour report used by the quadrature.
    pub fn contour(&self) -> &AbelJacobiContourReport {
        &self.contour
    }

    /// Returns the approximate integral value before reduction modulo `Λ`.
    pub fn value(&self) -> &Complex64 {
        &self.value
    }

    /// Returns the structured decomposition of the contour integral used by
    /// the current algorithm.
    pub fn decomposition(&self) -> &AbelJacobiIntegralDecomposition {
        &self.decomposition
    }

    /// Returns the initial square-root branch selected before continuity
    /// tracking along the contour begins.
    pub fn initial_branch_choice(&self) -> AbelJacobiInitialBranchChoice {
        self.decomposition.initial_branch_choice()
    }

    /// Returns the contribution of the finite segment from `X` to the anchor.
    pub fn segment_integral(&self) -> &Complex64 {
        self.decomposition.segment_integral()
    }

    /// Returns the contribution of the compactified outgoing ray before the
    /// asymptotic tail correction.
    pub fn ray_integral(&self) -> &Complex64 {
        self.decomposition.ray_integral()
    }

    /// Returns the asymptotic tail correction added beyond the sampled ray.
    pub fn tail_correction(&self) -> &Complex64 {
        self.decomposition.tail_correction()
    }

    /// Returns the quadrature step count used by the current run.
    pub fn integration_steps_used(&self) -> usize {
        self.numerics.integration_steps_used()
    }

    /// Returns how many branch-continuation corrections were applied.
    pub fn branch_adjustments_used(&self) -> usize {
        self.numerics.branch_adjustments_used()
    }

    /// Returns the comparison tolerance used by the current run.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.numerics.tolerance()
    }
}

/// Report for recovering one torus class from one point on an analytic
/// Weierstrass curve through an Abel-Jacobi integral.
///
/// The inverse map lands in the quotient `C / Λ`. Thus, the report keeps:
/// - the raw integral value before quotient reduction
/// - the canonical torus class after reduction modulo the recovered lattice
/// - one reduced representative in the chosen half-open fundamental
///   parallelogram
#[derive(Clone, Debug)]
pub struct AbelJacobiPointRecoveryReport {
    curve: AnalyticWeierstrassCurve,
    point: AnalyticCurvePoint,
    periods: RecoveredPeriodBasis,
    contour: AbelJacobiContourReport,
    raw_integral_value: Complex64,
    torus_point: ComplexTorusPoint,
    reduced_representative: Complex64,
    validation_report: AbelJacobiRoundtripValidationReport,
    metadata: AbelJacobiRecoveryMetadata,
}

impl AbelJacobiPointRecoveryReport {
    /// Builds one explicit Abel-Jacobi point-recovery report.
    pub(crate) fn new(
        periods: RecoveredPeriodBasis,
        integral_approx: AbelJacobiIntegralApprox,
        torus_point: ComplexTorusPoint,
        reduced_representative: Complex64,
        validation_report: AbelJacobiRoundtripValidationReport,
        metadata: AbelJacobiRecoveryMetadata,
    ) -> Result<Self, AnalyticCurveError> {
        let canonical_representative = periods
            .lattice()
            .point_from_fundamental_coordinates(torus_point.coordinate().clone());
        let comparison = ComplexApproxComparison::new(
            reduced_representative,
            canonical_representative,
            metadata.tolerance(),
        );

        if !comparison.agrees_approximately() {
            return Err(AnalyticCurveError::InverseUniformizationFailed);
        }

        Ok(Self {
            curve: integral_approx.curve().clone(),
            point: integral_approx.point().clone(),
            periods,
            contour: integral_approx.contour().clone(),
            raw_integral_value: *integral_approx.value(),
            torus_point,
            reduced_representative: canonical_representative,
            validation_report,
            metadata,
        })
    }

    /// Returns the source analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the source curve point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        &self.point
    }

    /// Returns the recovered period basis used to interpret the torus class.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        &self.periods
    }

    /// Returns the contour report used by the Abel-Jacobi quadrature.
    pub fn contour(&self) -> &AbelJacobiContourReport {
        &self.contour
    }

    /// Returns the implied period ratio `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.periods.tau()
    }

    /// Returns the raw Abel-Jacobi integral value before reduction modulo `Λ`.
    pub fn raw_integral_value(&self) -> &Complex64 {
        &self.raw_integral_value
    }

    /// Returns the recovered torus class as a canonical point of `C / Λ`.
    pub fn torus_point(&self) -> &ComplexTorusPoint {
        &self.torus_point
    }

    /// Returns one reduced complex representative of the recovered torus
    /// class.
    pub fn reduced_representative(&self) -> &Complex64 {
        &self.reduced_representative
    }

    /// Returns the forward-validation report attached to the recovered torus
    /// representative.
    pub fn validation_report(&self) -> &AbelJacobiRoundtripValidationReport {
        &self.validation_report
    }

    /// Returns the numerical metadata for the current run.
    pub fn metadata(&self) -> &AbelJacobiRecoveryMetadata {
        &self.metadata
    }
}

/// End-to-end report for point-level inverse uniformization from a curve-side
/// point `(x, y)` back to a torus class.
///
/// This report intentionally keeps the full period-recovery story visible
/// instead of exposing a second tau-only or lattice-only path:
///
/// 1. recover one period basis from the curve
/// 2. use that basis in the Abel-Jacobi inverse map
/// 3. package the resulting torus class together with its curve-side context
#[derive(Clone, Debug)]
pub struct InverseUniformizationPointRecoveryReport {
    period_basis_report: PeriodBasisRecoveryReport,
    point_recovery_report: AbelJacobiPointRecoveryReport,
}

impl InverseUniformizationPointRecoveryReport {
    /// Builds one explicit end-to-end inverse-uniformization report.
    pub(crate) fn new(
        period_basis_report: PeriodBasisRecoveryReport,
        point_recovery_report: AbelJacobiPointRecoveryReport,
    ) -> Result<Self, AnalyticCurveError> {
        if period_basis_report.curve() != point_recovery_report.curve() {
            return Err(AnalyticCurveError::InverseUniformizationFailed);
        }

        if period_basis_report.periods() != point_recovery_report.periods() {
            return Err(AnalyticCurveError::InverseUniformizationFailed);
        }

        Ok(Self {
            period_basis_report,
            point_recovery_report,
        })
    }

    /// Returns the full curve-level period-basis recovery report.
    pub fn period_basis_report(&self) -> &PeriodBasisRecoveryReport {
        &self.period_basis_report
    }

    /// Returns the Abel-Jacobi point-recovery layer.
    pub fn point_recovery_report(&self) -> &AbelJacobiPointRecoveryReport {
        &self.point_recovery_report
    }

    /// Returns the source analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        self.period_basis_report.curve()
    }

    /// Returns the source curve point.
    pub fn point(&self) -> &AnalyticCurvePoint {
        self.point_recovery_report.point()
    }

    /// Returns the recovered period basis.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.period_basis_report.periods()
    }

    /// Returns the implied period ratio `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.point_recovery_report.tau()
    }

    /// Returns the contour report used by the Abel-Jacobi quadrature.
    pub fn contour(&self) -> &AbelJacobiContourReport {
        self.point_recovery_report.contour()
    }

    /// Returns the recovered torus class.
    pub fn torus_point(&self) -> &ComplexTorusPoint {
        self.point_recovery_report.torus_point()
    }

    /// Returns one reduced complex representative of the recovered torus
    /// class.
    pub fn reduced_representative(&self) -> &Complex64 {
        self.point_recovery_report.reduced_representative()
    }

    /// Returns the forward-validation report attached to the recovered torus
    /// representative.
    pub fn validation_report(&self) -> &AbelJacobiRoundtripValidationReport {
        self.point_recovery_report.validation_report()
    }

    /// Returns the Abel-Jacobi numerical metadata.
    pub fn metadata(&self) -> &AbelJacobiRecoveryMetadata {
        self.point_recovery_report.metadata()
    }
}
