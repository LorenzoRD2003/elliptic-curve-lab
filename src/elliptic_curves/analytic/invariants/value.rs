use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, LatticeSumTruncation};
use crate::fields::complex_approx::ComplexApprox;

/// Approximate classical analytic invariants attached to a complex lattice `Λ`.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyticInvariants {
    g2: Complex64,
    g3: Complex64,
    discriminant: Complex64,
    j_invariant: Complex64,
    truncation: LatticeSumTruncation,
}

impl AnalyticInvariants {
    pub(crate) fn new(
        g2: Complex64,
        g3: Complex64,
        truncation: LatticeSumTruncation,
    ) -> Result<Self, AnalyticCurveError> {
        let discriminant = Self::discriminant_from_g2_g3(&g2, &g3);
        let j_invariant = Self::j_invariant_from_g2_g3(&g2, &g3)?;

        Ok(Self {
            g2,
            g3,
            discriminant,
            j_invariant,
            truncation,
        })
    }

    /// Approximation to `g₂(Λ) = 60 G₄(Λ)`.
    pub fn g2(&self) -> &Complex64 {
        &self.g2
    }

    /// Approximation to `g₃(Λ) = 140 G₆(Λ)`.
    pub fn g3(&self) -> &Complex64 {
        &self.g3
    }

    /// Approximation to `Δ(Λ) = g₂(Λ)^3 - 27 g₃(Λ)^2`.
    pub fn discriminant(&self) -> &Complex64 {
        &self.discriminant
    }

    /// Approximation to `j(Λ) = 1728 g₂(Λ)^3 / Δ(Λ)`.
    pub fn j_invariant(&self) -> &Complex64 {
        &self.j_invariant
    }

    /// Truncation policy used to compute `g₂` and `g₃`.
    pub fn truncation(&self) -> LatticeSumTruncation {
        self.truncation
    }

    /// Computes the classical discriminant expression `Δ = g₂^3 - 27 g₃^2`.
    pub fn discriminant_from_g2_g3(g2: &Complex64, g3: &Complex64) -> Complex64 {
        g2.powu(3) - Complex64::new(27.0, 0.0) * g3.powu(2)
    }

    /// Computes the classical analytic `j`-invariant from `g₂` and `g₃`.
    pub fn j_invariant_from_g2_g3(
        g2: &Complex64,
        g3: &Complex64,
    ) -> Result<Complex64, AnalyticCurveError> {
        let discriminant = Self::discriminant_from_g2_g3(g2, g3);

        if ComplexApprox::is_zero_with_tolerance(&discriminant, ComplexApprox::default_tolerance())
        {
            return Err(AnalyticCurveError::NearlySingularAnalyticCurve);
        }

        Ok(Complex64::new(1728.0, 0.0) * g2.powu(3) / discriminant)
    }
}
