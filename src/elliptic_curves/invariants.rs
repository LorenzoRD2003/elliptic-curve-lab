use crate::elliptic_curves::CurveModel;

/// Curve models that can expose a mathematically meaningful `j`-invariant.
///
/// This capability is intentionally narrower than [`CurveModel`]. Not every
/// curve presentation needs to expose classical invariants, and the crate
/// prefers small capability traits over inflating the base model trait.
///
/// The returned value should be invariant under isomorphism over an algebraic
/// closure in the sense appropriate to the represented curve family.
pub trait HasJInvariant: CurveModel {
    /// Returns the curve's `j`-invariant.
    fn j_invariant(&self) -> Self::Elem;
}
