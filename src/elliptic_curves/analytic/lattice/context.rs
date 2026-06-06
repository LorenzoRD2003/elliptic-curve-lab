use crate::elliptic_curves::analytic::{ComplexLattice, UpperHalfPlanePoint};

/// Shared access to the ambient upper-half-plane parameter and its standard
/// lattice.
#[allow(dead_code)]
pub(crate) trait HasAnalyticLatticeContext {
    /// Returns the upper-half-plane parameter `τ`.
    fn tau(&self) -> &UpperHalfPlanePoint;

    /// Returns the associated lattice `Λ_τ = ℤ + ℤτ`.
    fn lattice(&self) -> &ComplexLattice;
}
