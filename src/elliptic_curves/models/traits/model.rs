use crate::fields::traits::*;

/// Associates a base field, coordinate type, and point representation with a
/// curve model.
pub trait CurveModel {
    type Elem: Clone + core::fmt::Debug;
    type BaseField: Field<Elem = Self::Elem>;
    type Point: Clone + core::fmt::Debug;

    /// Returns whether the given point is the distinguished identity element.
    fn is_identity(&self, point: &Self::Point) -> bool;

    /// Returns whether the given point belongs to the curve model.
    fn contains(&self, point: &Self::Point) -> bool;

    /// Returns whether the point is on the curve and is not the identity.
    ///
    /// This is a small semantic convenience for APIs that want to talk about
    /// “finite curve points” or “non-identity points” without repeating the
    /// identity check at each call site.
    fn is_on_curve_nonzero(&self, point: &Self::Point) -> bool {
        self.contains(point) && !self.is_identity(point)
    }

    /// Returns the distinguished identity element of the curve model.
    ///
    /// For affine models this will typically be the explicit point at
    /// infinity. Other representations may encode the identity differently.
    fn identity(&self) -> Self::Point;
}

/// Minimal index sampler used by [`CurveModel::random_point`] without pulling
/// in an external randomness dependency.
///
/// The current crate intentionally avoids a `rand` dependency. This trait keeps
/// the curve-side API small while leaving room for callers to plug in a real
/// RNG adapter later if they want one.
pub trait PointIndexSampler {
    /// Chooses an index in `0..upper_bound`.
    ///
    /// Returning `None` lets the caller propagate sampling failure without
    /// introducing a dedicated randomness error type.
    fn sample_index(&mut self, upper_bound: usize) -> Option<usize>;
}

impl<T> PointIndexSampler for T
where
    T: FnMut(usize) -> Option<usize>,
{
    /// Delegates index selection to the wrapped closure.
    fn sample_index(&mut self, upper_bound: usize) -> Option<usize> {
        self(upper_bound)
    }
}
