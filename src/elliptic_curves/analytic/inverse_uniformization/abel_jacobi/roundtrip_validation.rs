use super::config::AbelJacobiConfig;
use super::report::AbelJacobiPointRecoveryReport;
use crate::ApproxTolerance;
use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::analytic::periods::RecoveredPeriodBasis;
use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, ComplexApproxComparison, EllipticFunctionTruncation,
    LatticeSumTruncation, map_torus_point_to_curve,
};

/// Explicit forward-validation policy for the point-level inverse-uniformization
/// roundtrip
///
/// `P -> z_P mod Λ -> (wp(z_P), wp'(z_P))`.
///
/// This is intentionally separate from the inverse Abel-Jacobi quadrature
/// budget. The point-recovery stage and the forward-validation stage are
/// different numerical problems, so callers should be able to tune them
/// independently.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointRoundTripValidationConfig {
    lattice_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    tolerance: ApproxTolerance,
}

impl PointRoundTripValidationConfig {
    /// Builds an explicit point-roundtrip validation policy.
    pub fn new(
        lattice_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        tolerance: ApproxTolerance,
    ) -> Self {
        Self {
            lattice_truncation,
            function_truncation,
            tolerance,
        }
    }

    /// Returns the baseline validation policy for educational experiments.
    pub fn educational_default() -> Self {
        Self::new(
            LatticeSumTruncation::new(16).expect("educational lattice truncation is valid"),
            EllipticFunctionTruncation::new(14).expect("educational function truncation is valid"),
            ApproxTolerance::educational_default(),
        )
    }

    /// Returns a tighter validation policy for more delicate runs.
    pub fn strict() -> Self {
        Self::new(
            LatticeSumTruncation::new(24).expect("strict lattice truncation is valid"),
            EllipticFunctionTruncation::new(22).expect("strict function truncation is valid"),
            ApproxTolerance::strict(),
        )
    }

    /// Returns a lighter validation policy for coarse exploratory work.
    pub fn loose() -> Self {
        Self::new(
            LatticeSumTruncation::new(12).expect("loose lattice truncation is valid"),
            EllipticFunctionTruncation::new(10).expect("loose function truncation is valid"),
            ApproxTolerance::loose(),
        )
    }

    /// Returns the lattice-sum truncation used by the forward validation.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.lattice_truncation
    }

    /// Returns the elliptic-function truncation used by the forward validation.
    pub fn function_truncation(&self) -> EllipticFunctionTruncation {
        self.function_truncation
    }

    /// Returns the tolerance used by the final `x`/`y` comparisons.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }
}

/// Successful forward-validation report for one recovered Abel-Jacobi torus
/// representative.
///
/// The report records how the recovered complex value behaved when it was sent
/// back through the direct torus-to-curve map `z -> (wp(z), wp'(z))`.
#[derive(Clone, Debug, PartialEq)]
pub struct AbelJacobiRoundtripValidationReport {
    recovered_curve_point: AnalyticCurvePoint,
    lattice_truncation: LatticeSumTruncation,
    function_truncation: EllipticFunctionTruncation,
    x_comparison: Option<ComplexApproxComparison>,
    y_comparison: Option<ComplexApproxComparison>,
}

impl AbelJacobiRoundtripValidationReport {
    /// Builds one explicit successful roundtrip-validation report.
    pub fn new(
        recovered_curve_point: AnalyticCurvePoint,
        lattice_truncation: LatticeSumTruncation,
        function_truncation: EllipticFunctionTruncation,
        x_comparison: Option<ComplexApproxComparison>,
        y_comparison: Option<ComplexApproxComparison>,
    ) -> Self {
        Self {
            recovered_curve_point,
            lattice_truncation,
            function_truncation,
            x_comparison,
            y_comparison,
        }
    }

    /// Returns the curve point recovered by the forward torus-to-curve map.
    pub fn recovered_curve_point(&self) -> &AnalyticCurvePoint {
        &self.recovered_curve_point
    }

    /// Returns the lattice-sum truncation used during the validation pass.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.lattice_truncation
    }

    /// Returns the elliptic-function truncation used during the validation pass.
    pub fn function_truncation(&self) -> EllipticFunctionTruncation {
        self.function_truncation
    }

    /// Returns the `x`-coordinate comparison when the validated point is
    /// finite on both sides.
    pub fn x_comparison(&self) -> Option<&ComplexApproxComparison> {
        self.x_comparison.as_ref()
    }

    /// Returns the `y`-coordinate comparison when the validated point is
    /// finite on both sides.
    pub fn y_comparison(&self) -> Option<&ComplexApproxComparison> {
        self.y_comparison.as_ref()
    }

    /// Returns the absolute residual in the recovered `x`-coordinate.
    ///
    /// For the point at infinity this is defined to be `0`.
    pub fn x_residual_norm(&self) -> f64 {
        self.x_comparison
            .as_ref()
            .map(ComplexApproxComparison::absolute_difference)
            .unwrap_or(0.0)
    }

    /// Returns the absolute residual in the recovered `y`-coordinate.
    ///
    /// For the point at infinity this is defined to be `0`.
    pub fn y_residual_norm(&self) -> f64 {
        self.y_comparison
            .as_ref()
            .map(ComplexApproxComparison::absolute_difference)
            .unwrap_or(0.0)
    }

    /// Returns whether the stored forward-validation comparisons agree at the
    /// requested tolerance.
    pub fn agrees_approximately(&self) -> bool {
        self.x_comparison
            .as_ref()
            .map(ComplexApproxComparison::agrees_approximately)
            .unwrap_or(true)
            && self
                .y_comparison
                .as_ref()
                .map(ComplexApproxComparison::agrees_approximately)
                .unwrap_or(true)
    }
}

/// Public end-to-end report for the point-level roundtrip experiment
///
/// `P -> z_P mod Λ -> (wp(z_P), wp'(z_P))`.
///
/// This wrapper keeps the full recovered torus-side data visible through the
/// embedded [`AbelJacobiPointRecoveryReport`], while also recording the
/// explicit forward-validation policy that was used for the final comparison
/// back on the curve.
#[derive(Clone, Debug)]
pub struct PointRoundTripValidationReport {
    point_recovery_report: AbelJacobiPointRecoveryReport,
    validation_config: PointRoundTripValidationConfig,
}

impl PointRoundTripValidationReport {
    /// Builds one explicit point-roundtrip validation report.
    pub fn new(
        point_recovery_report: AbelJacobiPointRecoveryReport,
        validation_config: PointRoundTripValidationConfig,
    ) -> Result<Self, AnalyticCurveError> {
        let validation_report = point_recovery_report.validation_report();

        if validation_report.lattice_truncation() != validation_config.lattice_truncation()
            || validation_report.function_truncation() != validation_config.function_truncation()
        {
            return Err(AnalyticCurveError::InverseUniformizationFailed);
        }

        for comparison in [
            validation_report.x_comparison(),
            validation_report.y_comparison(),
        ]
        .into_iter()
        .flatten()
        {
            if comparison.tolerance() != validation_config.tolerance() {
                return Err(AnalyticCurveError::InverseUniformizationFailed);
            }
        }

        Ok(Self {
            point_recovery_report,
            validation_config,
        })
    }

    /// Returns the embedded point-recovery layer `P -> z_P mod Λ`.
    pub fn point_recovery_report(&self) -> &AbelJacobiPointRecoveryReport {
        &self.point_recovery_report
    }

    /// Returns the final forward-validation layer
    /// `z_P mod Λ -> (wp(z_P), wp'(z_P))`.
    pub fn validation_report(&self) -> &AbelJacobiRoundtripValidationReport {
        self.point_recovery_report.validation_report()
    }

    /// Returns the source curve point `P`.
    pub fn point(&self) -> &AnalyticCurvePoint {
        self.point_recovery_report.point()
    }

    /// Returns the recovered torus class.
    pub fn torus_point(&self) -> &crate::elliptic_curves::analytic::ComplexTorusPoint {
        self.point_recovery_report.torus_point()
    }

    /// Returns one reduced complex representative of the recovered torus class.
    pub fn reduced_representative(&self) -> &num_complex::Complex64 {
        self.point_recovery_report.reduced_representative()
    }

    /// Returns the curve point obtained by mapping the recovered torus class
    /// back through `(wp, wp')`.
    pub fn recovered_curve_point(&self) -> &AnalyticCurvePoint {
        self.validation_report().recovered_curve_point()
    }

    /// Returns the `x`-coordinate comparison when both points are finite.
    pub fn x_comparison(&self) -> Option<&ComplexApproxComparison> {
        self.validation_report().x_comparison()
    }

    /// Returns the `y`-coordinate comparison when both points are finite.
    pub fn y_comparison(&self) -> Option<&ComplexApproxComparison> {
        self.validation_report().y_comparison()
    }

    /// Returns the absolute residual in the recovered `x`-coordinate.
    pub fn x_residual_norm(&self) -> f64 {
        self.validation_report().x_residual_norm()
    }

    /// Returns the absolute residual in the recovered `y`-coordinate.
    pub fn y_residual_norm(&self) -> f64 {
        self.validation_report().y_residual_norm()
    }

    /// Returns the forward-validation lattice-sum truncation.
    pub fn lattice_truncation(&self) -> LatticeSumTruncation {
        self.validation_config.lattice_truncation()
    }

    /// Returns the forward-validation elliptic-function truncation.
    pub fn function_truncation(&self) -> EllipticFunctionTruncation {
        self.validation_config.function_truncation()
    }

    /// Returns the tolerance used by the final `x`/`y` comparisons.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.validation_config.tolerance()
    }

    /// Returns whether the full point roundtrip agreed at the requested
    /// tolerance.
    pub fn agrees_approximately(&self) -> bool {
        self.validation_report().agrees_approximately()
    }
}

/// Validates one recovered Abel-Jacobi value by mapping it back through the
/// forward torus-to-curve uniformization and comparing against the source
/// curve point.
///
/// This is the main numerical sanity check for the inverse stage:
/// after recovering a complex value `z`, we evaluate the direct map
/// `z -> (wp(z), wp'(z))` with caller-configured forward-validation
/// truncations, then compare the recovered affine coordinates against the
/// original `(x, y)` within the requested tolerance.
///
/// The point at infinity is treated separately and must map back to the pole
/// case of the forward uniformization.
pub(super) fn point_roundtrip_validation_config_from_abel_config(
    config: AbelJacobiConfig,
) -> Result<PointRoundTripValidationConfig, AnalyticCurveError> {
    Ok(PointRoundTripValidationConfig::new(
        validation_lattice_truncation(config)?,
        validation_function_truncation(config)?,
        config.tolerance,
    ))
}

pub(super) fn point_roundtrip_validation_report_for_representative(
    point: &AnalyticCurvePoint,
    periods: &RecoveredPeriodBasis,
    representative: num_complex::Complex64,
    validation_config: PointRoundTripValidationConfig,
) -> Result<AbelJacobiRoundtripValidationReport, AnalyticCurveError> {
    let lattice_truncation = validation_config.lattice_truncation();
    let function_truncation = validation_config.function_truncation();
    let map_result = map_torus_point_to_curve(
        periods.lattice(),
        representative,
        lattice_truncation,
        function_truncation,
        validation_config.tolerance(),
    )?;

    match (point, map_result.point()) {
        (AffinePoint::Infinity, AffinePoint::Infinity) => {
            Ok(AbelJacobiRoundtripValidationReport::new(
                map_result.point().clone(),
                lattice_truncation,
                function_truncation,
                None,
                None,
            ))
        }
        (
            AffinePoint::Finite { x, y },
            AffinePoint::Finite {
                x: recovered_x,
                y: recovered_y,
            },
        ) => {
            let x_comparison =
                ComplexApproxComparison::new(*recovered_x, *x, validation_config.tolerance());
            let y_comparison =
                ComplexApproxComparison::new(*recovered_y, *y, validation_config.tolerance());

            Ok(AbelJacobiRoundtripValidationReport::new(
                map_result.point().clone(),
                lattice_truncation,
                function_truncation,
                Some(x_comparison),
                Some(y_comparison),
            ))
        }
        _ => Err(AnalyticCurveError::PeriodValidationFailed),
    }
}

/// Chooses the lattice-sum truncation used in the Abel-Jacobi validation
/// roundtrip.
///
/// This helper reads the explicit validation policy instead of deriving a
/// radius from the inverse quadrature budget.
fn validation_lattice_truncation(
    config: AbelJacobiConfig,
) -> Result<LatticeSumTruncation, AnalyticCurveError> {
    LatticeSumTruncation::new(config.validation_policy.lattice_truncation_radius)
}

/// Chooses the elliptic-function truncation used in the Abel-Jacobi
/// validation roundtrip.
///
/// As with [`validation_lattice_truncation`], this is controlled by the
/// explicit validation policy rather than being inferred from the inverse
/// quadrature budget.
fn validation_function_truncation(
    config: AbelJacobiConfig,
) -> Result<EllipticFunctionTruncation, AnalyticCurveError> {
    EllipticFunctionTruncation::new(config.validation_policy.function_truncation_radius)
}
