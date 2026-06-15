use crate::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve, LegendrePeriodIntegralReport, LegendreReduction,
    NumericalRecoveryMetadata, UpperHalfPlanePoint, WeierstrassCubicRoots,
    periods::RecoveredPeriodBasis,
};

/// A pedagogical report explaining a recovered period basis.
///
/// The intended story is:
///
/// - start from one explicit Legendre reduction
/// - compute the necessary complete elliptic integrals
/// - rescale those normalized periods back to the original curve
/// - package the resulting validated basis
///
/// The numerical recovery algorithm itself is intentionally still deferred.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveredPeriodBasisReport {
    reduction: LegendreReduction,
    integral_report: LegendrePeriodIntegralReport,
    basis: RecoveredPeriodBasis,
}

impl RecoveredPeriodBasisReport {
    /// Builds one recovered-period-basis report from explicit ingredients.
    pub(crate) fn new(
        reduction: LegendreReduction,
        integral_report: LegendrePeriodIntegralReport,
        basis: RecoveredPeriodBasis,
    ) -> Self {
        Self {
            reduction,
            integral_report,
            basis,
        }
    }

    /// Returns the Legendre reduction used by the report.
    pub fn reduction(&self) -> &LegendreReduction {
        &self.reduction
    }

    /// Returns the Legendre-side integral report used by the report.
    pub fn integral_report(&self) -> &LegendrePeriodIntegralReport {
        &self.integral_report
    }

    /// Returns the differential rescaling factor used to transport normalized
    /// Legendre periods back to the original curve.
    pub fn invariant_differential_scale(&self) -> num_complex::Complex64 {
        self.reduction.invariant_differential_scale()
    }

    /// Returns the recovered period basis itself.
    pub fn basis(&self) -> &RecoveredPeriodBasis {
        &self.basis
    }

    /// Returns the recovered period ratio `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.basis.tau()
    }
}

/// Complete curve-level report for one successful period-basis recovery run.
///
/// This bundles the main stages of the current recovery pipeline:
///
/// 1. recover the Weierstrass cubic roots from `g₂, g₃`
/// 2. choose one deterministic Legendre reduction
/// 3. evaluate `K(λ)` and `K(1 - λ)` through the AGM
/// 4. transport the normalized Legendre periods back to the original curve
///
/// This top-level report intentionally reuses the lower-level
/// [`RecoveredPeriodBasisReport`] instead of duplicating the Legendre-side
/// data again. Convenience accessors still expose the most important derived
/// views directly.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodBasisRecoveryReport {
    curve: AnalyticWeierstrassCurve,
    roots: WeierstrassCubicRoots,
    basis_report: RecoveredPeriodBasisReport,
    metadata: NumericalRecoveryMetadata,
}

impl PeriodBasisRecoveryReport {
    /// Builds one complete period-basis recovery report from explicit pieces.
    pub(crate) fn new(
        curve: AnalyticWeierstrassCurve,
        roots: WeierstrassCubicRoots,
        basis_report: RecoveredPeriodBasisReport,
        metadata: NumericalRecoveryMetadata,
    ) -> Self {
        Self {
            curve,
            roots,
            basis_report,
            metadata,
        }
    }

    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        &self.basis_report
    }

    pub fn legendre_reduction(&self) -> &LegendreReduction {
        self.basis_report.reduction()
    }

    pub fn k_report(&self) -> &LegendrePeriodIntegralReport {
        self.basis_report.integral_report()
    }

    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.basis_report.basis()
    }

    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.basis_report.tau()
    }

    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        &self.metadata
    }
}
