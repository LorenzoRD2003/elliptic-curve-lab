use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, AnalyticWeierstrassCurve, ComplexLattice,
    LatticeSumTruncation, ModularQParameter, UpperHalfPlanePoint,
    lattice::HasAnalyticLatticeContext, weierstrass_model::AnalyticShortWeierstrassModel,
};
use crate::fields::complex_approx::ComplexApprox;

/// Distinguishes the classical special upper-half-plane parameters used in the
/// internal analytic examples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SpecialTauKind {
    TauI,
    TauRho,
    Generic,
}

/// Distinguishes the classical special `j`-invariant regimes visible in the
/// first analytic experiments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SpecialJKind {
    Near1728,
    NearZero,
    Generic,
}

/// Aggregated internal report for one complex-analytic lattice experiment.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ComplexAnalyticCurveLabReport {
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
    pub(crate) fn from_tau(
        tau: UpperHalfPlanePoint,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let lattice = ComplexLattice::from_tau(tau.clone());
        let q_parameter = ModularQParameter::from_tau(tau.clone());
        let invariants = lattice.analytic_invariants(truncation)?;
        let analytic_curve = AnalyticWeierstrassCurve::new(*invariants.g2(), *invariants.g3())?;
        let short_model = AnalyticShortWeierstrassModel::from_analytic_curve(&analytic_curve);
        let special_tau_kind = classify_special_tau(&tau);
        let special_j_kind = classify_special_j(invariants.j_invariant());

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

    pub(crate) fn invariants(&self) -> &AnalyticInvariants {
        &self.invariants
    }

    pub(crate) fn analytic_curve(&self) -> &AnalyticWeierstrassCurve {
        &self.analytic_curve
    }

    pub(crate) fn short_model(&self) -> &AnalyticShortWeierstrassModel {
        &self.short_model
    }

    pub(crate) fn special_tau_kind(&self) -> SpecialTauKind {
        self.special_tau_kind
    }

    pub(crate) fn special_j_kind(&self) -> SpecialJKind {
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
