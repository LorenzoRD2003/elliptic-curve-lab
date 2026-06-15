use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticCurvePoint, AnalyticWeierstrassCurve, PeriodRecoveryConfig,
    RecoveredPeriodBasis,
    inverse_uniformization::abel_jacobi::{
        AbelJacobiConfig, AbelJacobiIntegralApprox, AbelJacobiPointRecoveryReport,
        InverseUniformizationPointRecoveryReport, PointRoundTripValidationConfig,
        PointRoundTripValidationReport, approximate_abel_jacobi_integral_impl,
        recover_torus_point_from_curve_point_impl,
        recover_torus_point_from_curve_point_with_periods_impl,
        validate_point_inverse_uniformization_roundtrip_impl,
        validate_point_inverse_uniformization_roundtrip_with_periods_impl,
    },
};

impl AnalyticWeierstrassCurve {
    /// Approximates the Abel-Jacobi integral attached to one point of this
    /// analytic Weierstrass curve.
    pub fn approximate_abel_jacobi_integral(
        &self,
        point: &AnalyticCurvePoint,
        config: AbelJacobiConfig,
    ) -> Result<AbelJacobiIntegralApprox, AnalyticCurveError> {
        approximate_abel_jacobi_integral_impl(self, point, config)
    }

    /// Recovers the torus class of one point on this analytic Weierstrass
    /// curve, assuming a period basis is already known.
    pub fn recover_torus_point_from_curve_point_with_periods(
        &self,
        point: &AnalyticCurvePoint,
        periods: &RecoveredPeriodBasis,
        config: AbelJacobiConfig,
    ) -> Result<AbelJacobiPointRecoveryReport, AnalyticCurveError> {
        recover_torus_point_from_curve_point_with_periods_impl(self, point, periods, config)
    }

    /// Recovers the torus class of one point on this analytic Weierstrass
    /// curve by first recovering a period basis and then applying the
    /// Abel-Jacobi inverse map.
    pub fn recover_torus_point_from_curve_point(
        &self,
        point: &AnalyticCurvePoint,
        period_config: PeriodRecoveryConfig,
        abel_jacobi_config: AbelJacobiConfig,
    ) -> Result<InverseUniformizationPointRecoveryReport, AnalyticCurveError> {
        recover_torus_point_from_curve_point_impl(self, point, period_config, abel_jacobi_config)
    }

    /// Validates one point-level inverse-uniformization roundtrip when a
    /// period basis is already known.
    pub fn validate_point_inverse_uniformization_roundtrip_with_periods(
        &self,
        point: &AnalyticCurvePoint,
        periods: &RecoveredPeriodBasis,
        abel_config: AbelJacobiConfig,
        validation_config: PointRoundTripValidationConfig,
    ) -> Result<PointRoundTripValidationReport, AnalyticCurveError> {
        validate_point_inverse_uniformization_roundtrip_with_periods_impl(
            self,
            point,
            periods,
            abel_config,
            validation_config,
        )
    }

    /// Validates one point-level inverse-uniformization roundtrip by first
    /// recovering a period basis from this curve.
    pub fn validate_point_inverse_uniformization_roundtrip(
        &self,
        point: &AnalyticCurvePoint,
        period_config: PeriodRecoveryConfig,
        abel_config: AbelJacobiConfig,
        validation_config: PointRoundTripValidationConfig,
    ) -> Result<PointRoundTripValidationReport, AnalyticCurveError> {
        validate_point_inverse_uniformization_roundtrip_impl(
            self,
            point,
            period_config,
            abel_config,
            validation_config,
        )
    }
}
