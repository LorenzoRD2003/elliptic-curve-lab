use num_bigint::{BigInt, BigUint};
use num_traits::{CheckedSub, One, Zero};

use crate::{
    elliptic_curves::{
        frobenius::cm::CmTraceCandidateError,
        traits::{BigScalarGroupCurveModel, EnumerableCurveModel, PointIndexSampler},
    },
    fields::traits::EnumerableFiniteField,
};

/// Curve models that can determine a CM trace sign from a point witness.
///
/// For a curve over `F_p`, a Frobenius trace `t` gives
/// `#E(F_p) = p + 1 - t`. Given an absolute candidate `|t|`, the method tests
/// the two possible group orders on a caller-supplied point `P`:
///
/// - if `[p + 1 - |t|]P = O` and `[p + 1 + |t|]P ≠ O`, it returns `+|t|`;
/// - if `[p + 1 + |t|]P = O` and `[p + 1 - |t|]P ≠ O`, it returns `-|t|`;
/// - if the point does not distinguish the two signs, it returns `None`.
///
/// The `None` case is intentional: a point of small order may be killed by both
/// candidate orders, and an unsuitable point may be killed by neither.
///
/// Complexity: `Θ(log(p + |t|))` group additions/doublings.
pub trait CmTraceSignCurveModel: BigScalarGroupCurveModel {
    /// Returns whether `[scalar]point = O` for the current curve model.
    ///
    /// This helper is intentionally scoped to the CM trace-sign story. It uses
    /// the same big-scalar multiplication surface as the sign test, rather
    /// than any small-scalar convenience on [`GroupCurveModel`].
    fn cm_scalar_kills_point(
        &self,
        point: &Self::Point,
        scalar: &BigUint,
    ) -> Result<bool, CmTraceCandidateError> {
        Ok(self
            .mul_scalar_biguint(point, scalar)
            .map(|multiple| self.is_identity(&multiple))?)
    }

    /// Determines the sign of one CM absolute-trace candidate using `point`.
    fn cm_trace_from_absolute_trace_with_point(
        &self,
        p: &BigUint,
        absolute_trace: &BigUint,
        point: &Self::Point,
    ) -> Result<Option<BigInt>, CmTraceCandidateError> {
        let base = p + BigUint::one();
        let Some(order_for_positive_trace) = base.checked_sub(absolute_trace) else {
            return Ok(None);
        };
        let order_for_negative_trace = base + absolute_trace;

        let positive_trace_order_kills_point =
            self.cm_scalar_kills_point(point, &order_for_positive_trace)?;
        let negative_trace_order_kills_point =
            self.cm_scalar_kills_point(point, &order_for_negative_trace)?;

        Ok(resolve_trace_sign_from_annihilation_tests(
            absolute_trace,
            positive_trace_order_kills_point,
            negative_trace_order_kills_point,
        ))
    }

    /// Enumerates rational points and returns the first witness that determines
    /// the sign of `absolute_trace`.
    ///
    /// This deterministic convenience is meant for the same small finite-field
    /// regime as [`EnumerableCurveModel`]. It tries non-identity finite points
    /// in enumeration order and returns `None` if none of them distinguishes
    /// the two candidate orders.
    ///
    /// Complexity: `Θ(n log(p + |t|))` group additions/doublings in the worst
    /// case, where `n` is the number of enumerated non-identity points.
    fn cm_trace_from_absolute_trace_by_enumeration(
        &self,
        p: &BigUint,
        absolute_trace: &BigUint,
    ) -> Result<Option<BigInt>, CmTraceCandidateError>
    where
        Self: EnumerableCurveModel,
        Self::BaseField: EnumerableFiniteField<Elem = Self::Elem>,
        Self::Point: PartialEq,
    {
        for point in self.finite_points() {
            if let Some(trace) =
                self.cm_trace_from_absolute_trace_with_point(p, absolute_trace, &point)?
            {
                return Ok(Some(trace));
            }
        }

        Ok(None)
    }

    /// Samples points and returns the first witness that determines the sign of
    /// `absolute_trace`.
    ///
    /// The crate intentionally avoids a direct `rand` dependency, so callers
    /// supply a [`PointIndexSampler`]. Each attempt samples one enumerated point
    /// through [`EnumerableCurveModel::random_point`] and applies
    /// [`Self::cm_trace_from_absolute_trace_with_point`]. Sampler exhaustion or
    /// `max_attempts` unsuccessful witnesses both return `None`.
    ///
    /// Complexity: `Θ(s log(p + |t|))` group additions/doublings plus `s`
    /// point-enumeration samples, where `s ≤ max_attempts`.
    fn cm_trace_from_absolute_trace_by_random_points<S>(
        &self,
        p: &BigUint,
        absolute_trace: &BigUint,
        sampler: &mut S,
        max_attempts: usize,
    ) -> Result<Option<BigInt>, CmTraceCandidateError>
    where
        Self: EnumerableCurveModel,
        Self::BaseField: EnumerableFiniteField<Elem = Self::Elem>,
        Self::Point: PartialEq,
        S: PointIndexSampler,
    {
        for _ in 0..max_attempts {
            let Some(point) = self.random_point(sampler) else {
                return Ok(None);
            };
            if let Some(trace) =
                self.cm_trace_from_absolute_trace_with_point(p, absolute_trace, &point)?
            {
                return Ok(Some(trace));
            }
        }

        Ok(None)
    }
}

impl<C> CmTraceSignCurveModel for C where C: BigScalarGroupCurveModel {}

fn resolve_trace_sign_from_annihilation_tests(
    absolute_trace: &BigUint,
    positive_trace_order_kills_point: bool,
    negative_trace_order_kills_point: bool,
) -> Option<BigInt> {
    match (
        positive_trace_order_kills_point,
        negative_trace_order_kills_point,
    ) {
        (true, false) => Some(BigInt::from(absolute_trace.clone())),
        (false, true) => Some(-BigInt::from(absolute_trace.clone())),
        (true, true) if absolute_trace.is_zero() => Some(BigInt::from(0u8)),
        _ => None,
    }
}
