use crate::fields::traits::FiniteField;

/// Finite fields whose full element set is small enough to enumerate directly.
///
/// This trait is intentionally narrower than [`FiniteField`]. It exists for
/// educational tasks such as exhaustive tables or curve-point enumeration, not
/// as a claim that every finite field backend should always materialize all of
/// its elements.
pub trait EnumerableFiniteField: FiniteField {
    /// Returns every field element in a deterministic order.
    fn elements() -> Vec<Self::Elem>;
}
