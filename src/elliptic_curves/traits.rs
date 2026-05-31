use crate::elliptic_curves::CurveError;
use crate::fields::{EnumerableFiniteField, Field, SqrtField};

/// Associates a base field, coordinate type, and point representation with a
/// curve model.
pub trait CurveModel {
    type Elem: Clone + core::fmt::Debug;
    type BaseField: Field<Elem = Self::Elem>;
    type Point;

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

/// Curve models that admit affine coordinate validation.
pub trait AffineCurveModel: CurveModel {
    /// Builds a point from affine coordinates after checking that it lies on
    /// the curve.
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError>;
}

/// Curve models equipped with an explicit additive group law on their points.
///
/// This trait is narrower than [`CurveModel`]: it is for concrete models where
/// the crate is ready to expose actual point addition, doubling, and scalar
/// multiplication, not just curve membership and constructors.
///
/// Implementations should treat `CurveError::PointNotOnCurve` as the honest
/// failure mode for finite inputs that do not belong to the model.
pub trait GroupCurveModel: CurveModel
where
    Self::Point: Clone,
{
    /// Returns the additive inverse of a point.
    ///
    /// For on-curve inputs this should stay on the curve and satisfy
    /// `P + (-P) = O`.
    fn neg(&self, point: &Self::Point) -> Self::Point;

    /// Adds two curve points.
    ///
    /// Implementations should return the honest group-law sum when both inputs
    /// lie on the curve, including the usual identity and inverse cases.
    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError>;

    /// Subtracts `right` from `left`.
    ///
    /// The default implementation negates `right` and then reuses
    /// [`Self::add`].
    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add(left, &negated)
    }

    /// Doubles a point under the group law.
    ///
    /// The default implementation calls [`Self::add`] with the same point
    /// twice after validating membership.
    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.add(point, point)
    }

    /// Multiplies a point by a non-negative integer using double-and-add.
    ///
    /// This is the clear educational baseline rather than an optimized
    /// cryptographic ladder.
    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        let mut result = self.identity();
        let mut base = point.clone();
        let mut k = scalar;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add(&result, &base)?;
            }

            k >>= 1;

            if k > 0 {
                base = self.double(&base)?;
            }
        }

        Ok(result)
    }

    /// Multiplies a point by a signed integer.
    ///
    /// Negative scalars are handled by negating the point and reusing
    /// [`Self::mul_scalar`].
    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar(point, scalar as u64)
        }
    }

    /// Returns whether the point is killed by multiplication by `n`.
    ///
    /// This helper is meant for the usual positive-integer notion of
    /// `n`-torsion. To avoid the degenerate convention `[0]P = O`, this method
    /// returns `false` when `n == 0`.
    ///
    /// Invalid off-curve inputs are treated honestly and return `false`.
    fn is_torsion_point(&self, point: &Self::Point, n: u64) -> bool {
        if n == 0 {
            return false;
        }

        self.mul_scalar(point, n)
            .is_ok_and(|multiple| self.is_identity(&multiple))
    }
}

/// Curve models that can lift an `x`-coordinate into one or two affine points.
///
/// This trait is intentionally about the *shape* of the curve equation rather
/// than about enumeration or group-law operations. It models the common
/// situation where the curve equation determines `y^2` from a chosen `x`.
pub trait LiftXCoordinate: AffineCurveModel
where
    Self::BaseField: SqrtField<Elem = Self::Elem>,
{
    /// Returns the right-hand side of the curve equation as a function of `x`.
    fn rhs(&self, x: &Self::Elem) -> Self::Elem;

    /// Builds one point above the given `x` when a square root exists.
    ///
    /// Which square root is chosen is delegated to the base field's
    /// [`SqrtField`] implementation.
    fn point_from_x(&self, x: Self::Elem) -> Option<Self::Point> {
        let y = Self::BaseField::sqrt(&self.rhs(&x))?;
        self.point(x, y).ok()
    }

    /// Builds the two points above the given `x` when square roots exist.
    ///
    /// When the only root is `0`, both returned points are the same because
    /// the two square roots coincide.
    fn points_from_x(&self, x: Self::Elem) -> Option<(Self::Point, Self::Point)> {
        let (left_y, right_y) = Self::BaseField::sqrt_pair(&self.rhs(&x))?;
        let left = self.point(x.clone(), left_y).ok()?;
        let right = self.point(x, right_y).ok()?;
        Some((left, right))
    }
}

/// Curve models that can be exhaustively enumerated over small finite fields.
///
/// This trait is intentionally narrower than [`CurveModel`]. It is meant for
/// educational scenarios where:
///
/// - the base field is small enough to enumerate directly
/// - the curve can be reconstructed by scanning every `x` and lifting
///   `y`-coordinates
///
/// It should not be read as a promise that every curve model in the crate
/// ought to support exhaustive point materialization.
pub trait EnumerableCurveModel: LiftXCoordinate
where
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem>,
    Self::Point: PartialEq,
{
    /// Returns all finite non-identity points under direct enumeration.
    ///
    /// The current algorithm enumerates every `x` in the base field, lifts the
    /// corresponding points, and deduplicates the `y = 0` case.
    fn finite_points(&self) -> Vec<Self::Point> {
        let mut points = Vec::new();

        for x in Self::BaseField::elements() {
            if let Some((left, right)) = self.points_from_x(x) {
                points.push(left);
                if points.last().is_some_and(|last| *last != right) {
                    points.push(right);
                }
            }
        }

        points
    }

    /// Returns every curve point, including the identity, under direct
    /// enumeration.
    fn points(&self) -> Vec<Self::Point> {
        let mut points = Vec::with_capacity(self.finite_points().len() + 1);
        points.push(self.identity());
        points.extend(self.finite_points());
        points
    }

    /// Returns the total number of points under direct enumeration.
    fn order(&self) -> usize {
        self.points().len()
    }

    /// Chooses one point using a minimal index-sampling interface.
    ///
    /// This samples from the fully enumerated point set. It is therefore meant
    /// only for the same small educational settings as [`EnumerableCurveModel`]
    /// itself.
    fn random_point<R>(&self, rng: &mut R) -> Option<Self::Point>
    where
        R: PointIndexSampler,
    {
        let mut points = self.points();
        let index = rng.sample_index(points.len())?;
        if index >= points.len() {
            return None;
        }

        Some(points.swap_remove(index))
    }
}

impl<T> EnumerableCurveModel for T
where
    T: LiftXCoordinate,
    T::BaseField: EnumerableFiniteField<Elem = T::Elem> + SqrtField<Elem = T::Elem>,
    T::Point: PartialEq,
{
}

/// Small finite curve groups whose full point set can be enumerated honestly.
///
/// This capability is narrower than both [`GroupCurveModel`] and
/// [`EnumerableCurveModel`]. It is meant for educational settings where the
/// whole curve group is small enough that group order and point order can be
/// computed directly.
///
/// In particular, helpers such as [`Self::point_order`] should be read as
/// small-group educational utilities backed by direct traversal, not as
/// efficient generic algorithms for large elliptic-curve groups.
pub trait FiniteGroupCurveModel: GroupCurveModel + EnumerableCurveModel
where
    Self::BaseField: EnumerableFiniteField<Elem = Self::Elem> + SqrtField<Elem = Self::Elem>,
    Self::Point: Clone + PartialEq,
{
    /// Returns the order of a point when the ambient curve group is small
    /// enough to enumerate directly.
    ///
    /// The order is the least positive integer `n` such that `[n]P = O`.
    /// Invalid off-curve inputs return `None` honestly. The identity has order
    /// `1`.
    ///
    /// The current implementation searches by repeated addition up to the full
    /// enumerated group order. That is intentional for small educational
    /// backends and should not be read as an efficient algorithm for large
    /// groups.
    fn point_order(&self, point: &Self::Point) -> Option<usize> {
        if !self.contains(point) {
            return None;
        }

        if self.is_identity(point) {
            return Some(1);
        }

        let mut accumulator = self.identity();
        for n in 1..=self.order() {
            accumulator = self.add(&accumulator, point).ok()?;
            if self.is_identity(&accumulator) {
                return Some(n);
            }
        }

        None
    }
}

impl<T> FiniteGroupCurveModel for T
where
    T: GroupCurveModel + EnumerableCurveModel,
    T::BaseField: EnumerableFiniteField<Elem = T::Elem> + SqrtField<Elem = T::Elem>,
    T::Point: Clone + PartialEq,
{
}
