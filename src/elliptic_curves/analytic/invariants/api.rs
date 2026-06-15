use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticCurveError, AnalyticInvariants, ComplexLattice, LatticeSumTruncation,
    UpperHalfPlanePoint,
};

impl ComplexLattice {
    /// Approximates the classical analytic invariant `g₂(Λ)`.
    pub fn analytic_g2(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<Complex64, AnalyticCurveError> {
        let g4 = self.g4_sum(truncation)?;
        Ok(Complex64::new(60.0, 0.0) * *g4.value())
    }

    /// Approximates the classical analytic invariant `g₃(Λ)`.
    pub fn analytic_g3(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<Complex64, AnalyticCurveError> {
        let g6 = self.g6_sum(truncation)?;
        Ok(Complex64::new(140.0, 0.0) * *g6.value())
    }

    /// Computes the approximate analytic invariants attached to `self`.
    pub fn analytic_invariants(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<AnalyticInvariants, AnalyticCurveError> {
        AnalyticInvariants::new(
            self.analytic_g2(truncation)?,
            self.analytic_g3(truncation)?,
            truncation,
        )
    }
}

impl UpperHalfPlanePoint {
    /// Computes the approximate analytic invariants attached to the standard
    /// lattice `Λ_τ = ℤ + ℤτ`.
    pub fn analytic_invariants(
        &self,
        truncation: LatticeSumTruncation,
    ) -> Result<AnalyticInvariants, AnalyticCurveError> {
        ComplexLattice::from_tau(self.clone()).analytic_invariants(truncation)
    }
}
