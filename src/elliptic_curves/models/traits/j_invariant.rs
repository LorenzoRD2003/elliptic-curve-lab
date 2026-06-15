use crate::elliptic_curves::traits::CurveModel;

/// Curve models that can expose a mathematically meaningful `j`-invariant.
///
/// The returned value should be invariant under isomorphism over an algebraic
/// closure in the sense appropriate to the represented curve family.
pub trait HasJInvariant: CurveModel {
    /// Returns the curve's `j`-invariant.
    fn j_invariant(&self) -> Self::Elem;
}
