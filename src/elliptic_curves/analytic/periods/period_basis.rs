use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticWeierstrassCurve, ComplexLattice, FundamentalDomainReductionReport,
    FundamentalDomainReductionStatus, UpperHalfPlanePoint,
    reduce_tau_to_standard_fundamental_domain,
};
use crate::elliptic_curves::analytic::periods::{
    LegendrePeriodIntegralReport, LegendreReduction, NumericalRecoveryMetadata,
    PeriodRecoveryConfig, PeriodRecoveryMethod, PeriodRecoveryStatus, WeierstrassCubicRoots,
    legendre_period_integral_report, recover_weierstrass_cubic_roots_with_report,
};
use crate::fields::ComplexApprox;

/// One recovered ordered period basis for an analytic elliptic curve.
///
/// This type intentionally wraps one validated [`ComplexLattice`] instead of
/// storing `ω₁`, `ω₂`, and `τ` as parallel fields. That keeps the non-
/// degeneracy and positive-orientation invariants in one place.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveredPeriodBasis {
    lattice: ComplexLattice,
}

impl RecoveredPeriodBasis {
    /// Builds a recovered period basis from explicit periods.
    pub fn new(omega1: Complex64, omega2: Complex64) -> Result<Self, AnalyticCurveError> {
        Ok(Self {
            lattice: ComplexLattice::new(omega1, omega2)?,
        })
    }

    /// Wraps an already validated complex lattice as a recovered period basis.
    pub fn from_lattice(lattice: ComplexLattice) -> Self {
        Self { lattice }
    }

    /// Returns the first recovered period `ω₁`.
    pub fn omega1(&self) -> &Complex64 {
        self.lattice.omega1()
    }

    /// Returns the second recovered period `ω₂`.
    pub fn omega2(&self) -> &Complex64 {
        self.lattice.omega2()
    }

    /// Returns the recovered period ratio `τ = ω₂ / ω₁`.
    ///
    /// Because this type only stores validated lattices, the ratio is expected
    /// to lie in the upper half-plane.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.lattice
            .tau()
            .expect("validated recovered period basis must have tau in the upper half-plane")
    }

    /// Returns the underlying validated lattice.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Consumes the wrapper and returns the underlying validated lattice.
    pub fn into_lattice(self) -> ComplexLattice {
        self.lattice
    }

    /// Returns the oriented area of the recovered period parallelogram.
    pub fn oriented_area(&self) -> f64 {
        self.lattice.oriented_area()
    }

    /// Returns the covolume of the recovered period lattice.
    pub fn covolume(&self) -> f64 {
        self.lattice.covolume()
    }
}

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
    pub fn new(
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
    pub fn invariant_differential_scale(&self) -> Complex64 {
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
    pub fn new(
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

    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the recovered Weierstrass cubic roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    /// Returns the lower-level Legendre-to-period-basis report.
    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        &self.basis_report
    }

    /// Returns the Legendre reduction used by the recovery pipeline.
    pub fn legendre_reduction(&self) -> &LegendreReduction {
        self.basis_report.reduction()
    }

    /// Returns the Legendre-side complete-elliptic-integral report.
    pub fn k_report(&self) -> &LegendrePeriodIntegralReport {
        self.basis_report.integral_report()
    }

    /// Returns the recovered period basis.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.basis_report.basis()
    }

    /// Returns the recovered period ratio `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.basis_report.tau()
    }

    /// Returns the aggregated numerical metadata for the whole recovery run.
    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        &self.metadata
    }
}

/// A focused curve-level report for users who primarily want the recovered
/// modular parameter `τ`.
///
/// This wrapper intentionally reuses [`PeriodBasisRecoveryReport`] instead of
/// re-running or re-encoding the recovery pipeline. It keeps the full period
/// basis available for callers who later decide they also need `ω₁`, `ω₂`, or
/// the Legendre-side diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct TauRecoveryReport {
    period_basis_report: PeriodBasisRecoveryReport,
}

impl TauRecoveryReport {
    /// Builds a `τ`-focused report from an already computed period-basis
    /// recovery report.
    pub fn new(period_basis_report: PeriodBasisRecoveryReport) -> Self {
        Self {
            period_basis_report,
        }
    }

    /// Returns the underlying complete period-basis recovery report.
    pub fn period_basis_report(&self) -> &PeriodBasisRecoveryReport {
        &self.period_basis_report
    }

    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        self.period_basis_report.curve()
    }

    /// Returns the recovered Weierstrass cubic roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        self.period_basis_report.roots()
    }

    /// Returns the lower-level Legendre-to-period-basis report.
    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        self.period_basis_report.basis_report()
    }

    /// Returns the Legendre reduction used by the recovery pipeline.
    pub fn legendre_reduction(&self) -> &LegendreReduction {
        self.period_basis_report.legendre_reduction()
    }

    /// Returns the Legendre-side complete-elliptic-integral report.
    pub fn k_report(&self) -> &LegendrePeriodIntegralReport {
        self.period_basis_report.k_report()
    }

    /// Returns the recovered period basis.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.period_basis_report.periods()
    }

    /// Returns the recovered modular parameter `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> UpperHalfPlanePoint {
        self.period_basis_report.tau()
    }

    /// Returns the aggregated numerical metadata for the whole recovery run.
    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        self.period_basis_report.metadata()
    }
}

/// A `τ`-recovery report together with an explicit canonicalization to the
/// standard fundamental domain of `SL₂(ℤ)`.
///
/// This report keeps two distinct layers visible:
///
/// - the naturally recovered `τ = ω₂ / ω₁` coming from one chosen period basis
/// - the canonical representative obtained by modular reduction
///
/// The second layer does not replace the first one. It explains how the
/// canonical representative is related to the naturally recovered value via an
/// accumulated modular matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct CanonicalTauRecoveryReport {
    tau_recovery_report: TauRecoveryReport,
    fundamental_domain_reduction: FundamentalDomainReductionReport,
}

impl CanonicalTauRecoveryReport {
    /// Builds a canonical-`τ` recovery report from its two explicit layers.
    pub fn new(
        tau_recovery_report: TauRecoveryReport,
        fundamental_domain_reduction: FundamentalDomainReductionReport,
    ) -> Self {
        Self {
            tau_recovery_report,
            fundamental_domain_reduction,
        }
    }

    /// Returns the underlying natural `τ`-recovery report.
    pub fn tau_recovery_report(&self) -> &TauRecoveryReport {
        &self.tau_recovery_report
    }

    /// Returns the modular reduction report that canonicalizes the recovered
    /// `τ`.
    pub fn fundamental_domain_reduction(&self) -> &FundamentalDomainReductionReport {
        &self.fundamental_domain_reduction
    }

    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        self.tau_recovery_report.curve()
    }

    /// Returns the recovered Weierstrass cubic roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        self.tau_recovery_report.roots()
    }

    /// Returns the lower-level Legendre-to-period-basis report.
    pub fn basis_report(&self) -> &RecoveredPeriodBasisReport {
        self.tau_recovery_report.basis_report()
    }

    /// Returns the recovered period basis before modular canonicalization.
    pub fn periods(&self) -> &RecoveredPeriodBasis {
        self.tau_recovery_report.periods()
    }

    /// Returns the naturally recovered modular parameter
    /// `τ = ω₂ / ω₁` before any `SL₂(ℤ)` normalization.
    pub fn original_tau(&self) -> UpperHalfPlanePoint {
        self.tau_recovery_report.tau()
    }

    /// Returns the canonical representative in the standard fundamental
    /// domain.
    pub fn canonical_tau(&self) -> &UpperHalfPlanePoint {
        self.fundamental_domain_reduction.reduced_tau()
    }

    /// Returns the accumulated modular matrix `γ` with
    /// `canonical_tau = γ(original_tau)`.
    pub fn accumulated_matrix(&self) -> crate::elliptic_curves::ModularMatrix {
        self.fundamental_domain_reduction.accumulated_matrix()
    }

    /// Returns the aggregated numerical metadata for the original recovery
    /// pipeline.
    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        self.tau_recovery_report.metadata()
    }
}

/// Recovers one period-basis report from an explicit Legendre reduction.
///
/// If the chosen Legendre normalization is
///
/// `Y² = X(X-1)(X-λ)`,
///
/// then the standard complete-elliptic-integral values
///
/// `2 K(λ)` and `2 i K(1-λ)`.
///
/// These are the classical **half-period** integrals for the Legendre model.
/// The full period lattice of `℘` is generated by twice those values.
///
/// The original Weierstrass differential satisfies
///
/// `dx / y = scale * dX / Y`,
///
/// where `scale = reduction.invariant_differential_scale()`.
///
/// So this routine recovers the concrete **full** period basis
///
/// `ω₁ = 4 scale K(λ)`,
/// `ω₂ = 4 scale i K(1-λ)`.
///
/// It then validates that the induced ratio `τ = ω₂ / ω₁` still agrees with
/// the already computed Legendre-side candidate and lies in the upper
/// half-plane.
///
/// Complexity: `Θ(a)`, where `a = config.agm_max_iterations()`.
/// The dominant work comes from the two AGM evaluations inside
/// [`legendre_period_integral_report`].
pub fn recover_period_basis_from_legendre_reduction(
    reduction: &LegendreReduction,
    config: PeriodRecoveryConfig,
) -> Result<RecoveredPeriodBasisReport, AnalyticCurveError> {
    let integral_report = legendre_period_integral_report(reduction, config)?;
    let basis = build_period_basis_from_legendre_data(reduction, &integral_report, config)?;

    Ok(RecoveredPeriodBasisReport::new(
        reduction.clone(),
        integral_report,
        basis,
    ))
}

/// Recovers one complete period-basis report directly from an analytic
/// Weierstrass curve.
///
/// 1. recover the cubic roots of `4x³ - g₂x - g₃`
/// 2. choose a deterministic Legendre representative away from `{0,1,∞}`
/// 3. evaluate `K(λ)` and `K(1 - λ)` via the complex AGM
/// 4. rescale the normalized Legendre periods back to the original curve
///
/// The resulting metadata aggregates the Newton work from cubic-root recovery
/// with the AGM work from both complete elliptic integrals.
///
/// Complexity: `Θ(n + a)`, where
/// - `n = config.newton_max_iterations()`
/// - `a = config.agm_max_iterations()`
pub fn recover_period_basis(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<PeriodBasisRecoveryReport, AnalyticCurveError> {
    let root_report = recover_weierstrass_cubic_roots_with_report(curve, config)?;
    let roots = root_report.roots().clone();
    let legendre_reduction = LegendreReduction::from_roots(&roots, config.tolerance())?;
    let basis_report = recover_period_basis_from_legendre_reduction(&legendre_reduction, config)?;
    let metadata = build_complete_period_recovery_metadata(
        root_report.metadata(),
        basis_report.integral_report(),
        *basis_report.tau().tau(),
    );

    Ok(PeriodBasisRecoveryReport::new(
        curve.clone(),
        roots,
        basis_report,
        metadata,
    ))
}

/// Recovers the modular parameter `τ` directly from an analytic Weierstrass
/// curve.
///
/// This is a convenience wrapper around [`recover_period_basis`]. It does not
/// introduce a second numerical route: the same cubic-root recovery, Legendre
/// reduction, AGM-based complete elliptic integrals, and period transport are
/// reused exactly once, then re-exposed through a `τ`-focused report.
///
/// Complexity: `Θ(n + a)`, where
/// - `n = config.newton_max_iterations()`
/// - `a = config.agm_max_iterations()`
pub fn recover_tau_from_curve(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<TauRecoveryReport, AnalyticCurveError> {
    recover_period_basis(curve, config).map(TauRecoveryReport::new)
}

/// Recovers a canonically normalized modular parameter from an analytic
/// Weierstrass curve.
///
/// This function deliberately composes two mathematically distinct stages.
///
/// 1. Recover one natural modular parameter `τ = ω₂ / ω₁` from the period
///    basis attached to the original curve.
/// 2. Reduce that recovered `τ` to the standard fundamental domain of
///    `SL₂(ℤ)`.
///
/// Keeping those stages separate is important pedagogically:
///
/// - [`recover_tau_from_curve`] exposes the `τ` value induced directly by the
///   recovered period basis.
/// - this function then explains how to choose a canonical modular
///   representative of the same lattice class.
///
/// The modular reduction uses
/// `config.fundamental_domain_reduction_max_steps()` transformations. If that
/// budget is hit before reaching the standard domain, this routine returns
/// [`AnalyticCurveError::PeriodValidationFailed`] instead of pretending that a
/// canonical representative was obtained.
///
/// Complexity: `Θ(n + a + m)`, where
///
/// - `n = config.newton_max_iterations()`
/// - `a = config.agm_max_iterations()`
/// - `m = config.fundamental_domain_reduction_max_steps()` is the modular
///   reduction step budget
///
/// In practice the recovery work dominates and the modular reduction stage is
/// tiny.
pub fn recover_canonical_tau_from_curve(
    curve: &AnalyticWeierstrassCurve,
    config: PeriodRecoveryConfig,
) -> Result<CanonicalTauRecoveryReport, AnalyticCurveError> {
    let tau_recovery_report = recover_tau_from_curve(curve, config)?;
    let fundamental_domain_reduction =
        reduce_recovered_tau_to_standard_domain(&tau_recovery_report, config)?;

    Ok(CanonicalTauRecoveryReport::new(
        tau_recovery_report,
        fundamental_domain_reduction,
    ))
}

/// Transports the normalized Legendre periods back to the original Weierstrass curve.
///
/// Suppose the chosen Legendre reduction writes the original cubic in the form
///
/// `x = e₂ + (e₁ - e₂) X`,
/// `y = legendre_y_scale * Y`,
/// `Y² = X(X-1)(X-λ)`.
///
/// For the invariant differential, the reduction stores the exact scale relation
/// `dx / y = differential_scale * dX / Y`,
/// where `differential_scale = reduction.invariant_differential_scale()`.
///
/// On the normalized Legendre side, our current AGM convention produces the
/// standard complete-elliptic-integral values
///
/// `2 K(λ)` and `2i K(1 -λ )`
///
/// for the differential `dX / Y`.
///
/// In the classical Weierstrass-to-Legendre bridge these are half-period
/// integrals: the full period lattice of `℘` is generated by twice those
/// values. Multiplying by `differential_scale` and doubling once more
/// therefore yields the candidate full periods for the original curve:
///
/// `ω₁ = 4 * differential_scale * K(λ)`,
/// `ω₂ = 4i * differential_scale * K(1 - λ)`.
///
/// This helper performs exactly that transport, then applies two numerical
/// sanity checks:
///
/// 1. both candidate periods must remain finite complex numbers
/// 2. the induced ratio `τ = ω₂ / ω₁` must agree, up to the configured
///    tolerance, with the already computed Legendre-side quantity
///    `tau_candidate = i K(1 - λ) / K(λ)`
///
/// The second check is mathematically natural: the common differential scale
/// cancels in the ratio, so the basis recovered on the original curve should
/// induce the same modular parameter as the normalized Legendre model.
///
/// Complexity: `Θ(1)` once the Legendre reduction and complete elliptic
/// integrals have already been computed.
fn build_period_basis_from_legendre_data(
    reduction: &LegendreReduction,
    integral_report: &LegendrePeriodIntegralReport,
    config: PeriodRecoveryConfig,
) -> Result<RecoveredPeriodBasis, AnalyticCurveError> {
    let differential_scale = reduction.invariant_differential_scale();
    let omega1 = Complex64::new(4.0, 0.0) * differential_scale * *integral_report.k_lambda.value();
    let omega2 =
        Complex64::new(0.0, 4.0) * differential_scale * *integral_report.k_complementary.value();

    if !omega1.is_finite() || !omega2.is_finite() {
        return Err(AnalyticCurveError::PeriodRecoveryFailed);
    }

    let basis = RecoveredPeriodBasis::new(omega1, omega2).map_err(map_basis_construction_error)?;
    validate_period_basis_against_tau_candidate(&basis, integral_report, config)?;

    Ok(basis)
}

fn validate_period_basis_against_tau_candidate(
    basis: &RecoveredPeriodBasis,
    integral_report: &LegendrePeriodIntegralReport,
    config: PeriodRecoveryConfig,
) -> Result<(), AnalyticCurveError> {
    let tau = basis.tau();
    let tau_value = *tau.tau();
    if !ComplexApprox::eq_with_tolerance(
        &tau_value,
        &integral_report.tau_candidate,
        config.tolerance(),
    ) {
        return Err(AnalyticCurveError::PeriodValidationFailed);
    }

    Ok(())
}

fn reduce_recovered_tau_to_standard_domain(
    tau_recovery_report: &TauRecoveryReport,
    config: PeriodRecoveryConfig,
) -> Result<FundamentalDomainReductionReport, AnalyticCurveError> {
    let reduction = reduce_tau_to_standard_fundamental_domain(
        tau_recovery_report.tau(),
        config.fundamental_domain_reduction_max_steps(),
    )?;

    match reduction.status() {
        FundamentalDomainReductionStatus::AlreadyReduced
        | FundamentalDomainReductionStatus::Reduced => Ok(reduction),
        FundamentalDomainReductionStatus::StepLimitReached => {
            Err(AnalyticCurveError::PeriodValidationFailed)
        }
    }
}

fn build_complete_period_recovery_metadata(
    root_metadata: &NumericalRecoveryMetadata,
    integral_report: &LegendrePeriodIntegralReport,
    recovered_tau: Complex64,
) -> NumericalRecoveryMetadata {
    let agm_iterations_used = integral_report.k_lambda.metadata().agm_iterations_used()
        + integral_report
            .k_complementary
            .metadata()
            .agm_iterations_used();
    let status = combine_period_basis_status(root_metadata, integral_report);
    let validation_residual_norm = Some(
        root_metadata
            .validation_residual_norm()
            .unwrap_or(0.0)
            .max((recovered_tau - integral_report.tau_candidate).norm()),
    );

    let mut metadata = NumericalRecoveryMetadata::new(
        PeriodRecoveryMethod::Hybrid,
        status,
        root_metadata.newton_iterations_used(),
        agm_iterations_used,
        0,
        root_metadata.branch_lattice_searches_used(),
        root_metadata.tolerance(),
        validation_residual_norm,
    );

    if let (
        Some(cardano_discriminant),
        Some(cardano_product_residual_norm),
        Some(selected_u_branch_index),
        Some(selected_v_branch_index),
    ) = (
        root_metadata.cardano_discriminant().copied(),
        root_metadata.cardano_product_residual_norm(),
        root_metadata.selected_u_branch_index(),
        root_metadata.selected_v_branch_index(),
    ) {
        metadata = metadata.with_cardano_diagnostics(
            cardano_discriminant,
            cardano_product_residual_norm,
            selected_u_branch_index,
            selected_v_branch_index,
        );
    }

    metadata
}

fn combine_period_basis_status(
    root_metadata: &NumericalRecoveryMetadata,
    integral_report: &LegendrePeriodIntegralReport,
) -> PeriodRecoveryStatus {
    [
        root_metadata.status(),
        integral_report.k_lambda.metadata().status(),
        integral_report.k_complementary.metadata().status(),
    ]
    .into_iter()
    .find(|status| *status != PeriodRecoveryStatus::Succeeded)
    .unwrap_or(PeriodRecoveryStatus::Succeeded)
}

fn map_basis_construction_error(error: AnalyticCurveError) -> AnalyticCurveError {
    match error {
        AnalyticCurveError::NonPositiveLatticeOrientation => {
            AnalyticCurveError::PeriodRatioNotInUpperHalfPlane
        }
        AnalyticCurveError::DegenerateLattice => AnalyticCurveError::PeriodRecoveryFailed,
        other => other,
    }
}
