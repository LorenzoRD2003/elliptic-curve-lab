use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, ComplexLattice, HasAnalyticLatticeContext, UpperHalfPlanePoint,
};

/// One chosen approximate period basis for the analytic uniformization
/// lattice of a complex elliptic curve.
///
/// The same curve admits many ordered bases related by `SL₂(ℤ)`, so this type
/// intentionally stores one chosen basis rather than pretending to represent a
/// canonical pair of periods. The cached modulus `τ = ω₂ / ω₁` is included
/// because it is the natural parameter for later modular-normalization and
/// `j`-comparison experiments.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodLatticeApprox {
    lattice: ComplexLattice,
    tau: UpperHalfPlanePoint,
}

impl PeriodLatticeApprox {
    /// Builds a period-lattice approximation from a validated lattice basis.
    pub fn new(lattice: ComplexLattice) -> Result<Self, AnalyticCurveError> {
        let tau = lattice.tau()?;
        Ok(Self { lattice, tau })
    }

    /// Builds the standard normalized basis `ω₁ = 1`, `ω₂ = τ`.
    ///
    /// This does not recover periods from a curve. It only packages the
    /// canonical lattice representative `Λ_τ = ℤ + ℤτ` when the modulus `τ`
    /// is already known.
    pub fn standard_from_tau(tau: UpperHalfPlanePoint) -> Self {
        Self {
            lattice: ComplexLattice::from_tau(tau.clone()),
            tau,
        }
    }

    /// Returns the stored lattice basis.
    pub fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }

    /// Returns the first period `ω₁`.
    pub fn omega1(&self) -> &Complex64 {
        self.lattice.omega1()
    }

    /// Returns the second period `ω₂`.
    pub fn omega2(&self) -> &Complex64 {
        self.lattice.omega2()
    }

    /// Returns the associated modulus `τ = ω₂ / ω₁`.
    pub fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }
}

impl HasAnalyticLatticeContext for PeriodLatticeApprox {
    fn tau(&self) -> &UpperHalfPlanePoint {
        &self.tau
    }

    fn lattice(&self) -> &ComplexLattice {
        &self.lattice
    }
}
