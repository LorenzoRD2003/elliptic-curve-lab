use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, AnalyticShortWeierstrassModel,
    AnalyticWeierstrassCurve, ComplexLattice, HasAnalyticLatticeContext, LatticeSumTruncation,
    ModularQParameter, UpperHalfPlanePoint, analytic_invariants,
};
use crate::fields::ComplexApprox;

/// Distinguishes the classical special upper-half-plane parameters used
/// throughout the analytic examples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpecialTauKind {
    /// The square-lattice point `τ = i`.
    TauI,
    /// The equianharmonic point `τ = ρ = -1/2 + (√3/2)i`.
    TauRho,
    /// Any other currently non-special sample point.
    Generic,
}

/// Distinguishes the classical special `j`-invariant regimes visible in the
/// first analytic experiments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpecialJKind {
    /// `j` is numerically close to `1728`, as for the square lattice.
    Near1728,
    /// `j` is numerically close to `0`, as for the equianharmonic lattice.
    NearZero,
    /// No special cancellation is currently detected.
    Generic,
}

/// Aggregated educational report for one complex-analytic lattice experiment.
///
/// This bundles together the main objects computed from a single parameter
/// `τ`:
///
/// - the standard lattice `Λ_τ = ℤ + ℤτ`
/// - the modular cusp parameter `q = e^{2π i τ}`
/// - the approximate analytic invariants `g₂`, `g₃`, `Δ`, and `j`
/// - the analytic cubic `y² = 4x³ - g₂x - g₃`
/// - its short-Weierstrass companion `y² = x³ + ax + b`
/// - lightweight classifications of the special `τ` and `j` cases
///
/// The point of this report is not to add new mathematics, but to keep the
/// first analytic layer coherent as one inspectable bundle.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexAnalyticCurveLabReport {
    tau: UpperHalfPlanePoint,
    lattice: ComplexLattice,
    q_parameter: ModularQParameter,
    invariants: AnalyticInvariants,
    analytic_curve: AnalyticWeierstrassCurve,
    short_model: AnalyticShortWeierstrassModel,
    special_tau_kind: SpecialTauKind,
    special_j_kind: SpecialJKind,
}

impl ComplexAnalyticCurveLabReport {
    /// Builds the aggregated analytic report attached to one `τ`.
    ///
    /// The only numerical knob here is the lattice truncation used to compute
    /// the approximate invariants; the `q`-parameter itself is exact at the
    /// current `Complex64` level.
    pub fn from_tau(
        tau: UpperHalfPlanePoint,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let lattice = ComplexLattice::from_tau(tau.clone());
        let q_parameter = ModularQParameter::from_tau(tau.clone());
        let invariants = analytic_invariants(&lattice, truncation)?;
        let analytic_curve = AnalyticWeierstrassCurve::new(invariants.g2, invariants.g3)?;
        let short_model = AnalyticShortWeierstrassModel::from_analytic_curve(&analytic_curve);
        let special_tau_kind = classify_special_tau(&tau);
        let special_j_kind = classify_special_j(&invariants.j_invariant);

        Ok(Self {
            tau,
            lattice,
            q_parameter,
            invariants,
            analytic_curve,
            short_model,
            special_tau_kind,
            special_j_kind,
        })
    }

    /// Returns the original upper-half-plane parameter `τ`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    /// Returns the standard lattice `Λ_τ = ℤ + ℤτ`.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the derived modular parameter `q = e^{2π i τ}`.
    pub fn q_parameter(&self) -> &ModularQParameter {
        &self.q_parameter
    }

    /// Returns the approximate analytic invariants.
    pub fn invariants(&self) -> &AnalyticInvariants {
        &self.invariants
    }

    /// Returns the analytic cubic `y² = 4x³ - g₂x - g₃`.
    pub fn analytic_curve(&self) -> &AnalyticWeierstrassCurve {
        &self.analytic_curve
    }

    /// Returns the short-Weierstrass companion model.
    pub fn short_model(&self) -> &AnalyticShortWeierstrassModel {
        &self.short_model
    }

    /// Returns the lightweight classification of the supplied `τ`.
    pub fn special_tau_kind(&self) -> SpecialTauKind {
        self.special_tau_kind
    }

    /// Returns the lightweight classification of the computed `j`.
    pub fn special_j_kind(&self) -> SpecialJKind {
        self.special_j_kind
    }
}

impl HasAnalyticLatticeContext for ComplexAnalyticCurveLabReport {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}

fn classify_special_tau(tau: &UpperHalfPlanePoint) -> SpecialTauKind {
    if ComplexApprox::eq_with_tolerance(
        tau.tau(),
        UpperHalfPlanePoint::tau_i().tau(),
        ComplexApprox::default_tolerance(),
    ) {
        return SpecialTauKind::TauI;
    }

    if ComplexApprox::eq_with_tolerance(
        tau.tau(),
        UpperHalfPlanePoint::tau_rho().tau(),
        ComplexApprox::default_tolerance(),
    ) {
        return SpecialTauKind::TauRho;
    }

    SpecialTauKind::Generic
}

fn classify_special_j(j_invariant: &Complex64) -> SpecialJKind {
    if ComplexApprox::eq_with_tolerance(
        j_invariant,
        &Complex64::new(1728.0, 0.0),
        crate::numerics::ApproxTolerance::new(1.0e-3, 1.0e-6),
    ) {
        return SpecialJKind::Near1728;
    }

    if ComplexApprox::eq_with_tolerance(
        j_invariant,
        &Complex64::new(0.0, 0.0),
        crate::numerics::ApproxTolerance::new(1.0e-4, 1.0e-6),
    ) {
        return SpecialJKind::NearZero;
    }

    SpecialJKind::Generic
}
