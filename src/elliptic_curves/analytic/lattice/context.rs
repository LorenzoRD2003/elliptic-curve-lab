use crate::elliptic_curves::analytic::{ComplexLattice, UpperHalfPlanePoint};

/// Shared access to the ambient upper-half-plane parameter and its standard
/// lattice.
pub trait HasAnalyticLatticeContext {
    /// Returns the upper-half-plane parameter `τ`.
    fn tau(&self) -> &UpperHalfPlanePoint;

    /// Returns the associated lattice `Λ_τ = ℤ + ℤτ`.
    fn lattice(&self) -> &ComplexLattice;
}
