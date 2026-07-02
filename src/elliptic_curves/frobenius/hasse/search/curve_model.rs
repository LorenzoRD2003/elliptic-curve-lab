use std::hash::Hash;

use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError,
    frobenius::{
        HasseInterval,
        hasse::search::{
            HasseBsgsConfig, HasseMultipleSearchReport,
            bsgs::find_annihilating_multiple_in_interval_bsgs_with_config_impl,
            naive::find_annihilating_multiple_in_interval_naive_report,
        },
    },
    traits::GroupCurveModel,
};

/// Curve models with an additive group law that can search one given Hasse
/// interval by either the naive scan or the BSGS refinement.
///
/// This trait is crate-private because it is internal execution machinery for
/// Frobenius/Hasse workflows rather than a stable user-facing abstraction.
pub(crate) trait HasseIntervalSearchCurveModel: GroupCurveModel
where
    Self::Point: Clone,
{
    /// Searches one already-chosen interval from left to right until
    /// `[M]P = O` is found or the interval is exhausted.
    ///
    /// Complexity: one `BigUint` scalar multiplication to build `[L]P`, then
    /// `Θ(|H|)` group additions, where `|H|` is the number of integer
    /// candidates in the supplied interval.
    fn find_annihilating_multiple_in_interval_naive(
        &self,
        point: &Self::Point,
        interval: HasseInterval,
    ) -> Result<HasseMultipleSearchReport<Self::Point>, CurveError> {
        find_annihilating_multiple_in_interval_naive_report(self, point, interval)
    }

    /// Searches one already-chosen interval with the baby-step/giant-step
    /// method from Algorithm 7.9 in the MIT 18.783 notes.
    ///
    /// This helper returns one `M ∈ H(q)` with `[M]P = O`, if found.
    ///
    /// Complexity: Let `c = |H(q) ∩ Z|`. The current implementation chooses
    /// `r = ceil(√c)` and `s = ceil(c/r)`, then performs:
    ///
    /// - `Θ(r)` group additions to build the baby steps
    /// - `Θ(1)` big-scalar multiplications to build `[a]P` and `[r]P`
    /// - `Θ(s)` hash lookups and giant-step additions
    ///
    /// Thus the dominant group-operation count is `Θ(r + s) = Θ(√c)`,
    /// which for Hasse intervals is `Θ(∜q)`.
    #[allow(dead_code)]
    fn find_annihilating_multiple_in_interval_bsgs(
        &self,
        point: &Self::Point,
        interval: HasseInterval,
    ) -> Result<Option<BigUint>, CurveError>
    where
        Self::Point: Eq + Hash,
    {
        self.find_annihilating_multiple_in_interval_bsgs_with_config(
            point,
            interval,
            HasseBsgsConfig::default(),
        )
    }

    /// Internal configurable BSGS engine for one Hasse interval.
    fn find_annihilating_multiple_in_interval_bsgs_with_config(
        &self,
        point: &Self::Point,
        interval: HasseInterval,
        config: HasseBsgsConfig,
    ) -> Result<Option<BigUint>, CurveError>
    where
        Self::Point: Eq + Hash,
    {
        find_annihilating_multiple_in_interval_bsgs_with_config_impl(self, point, interval, config)
    }
}

impl<C: GroupCurveModel + ?Sized> HasseIntervalSearchCurveModel for C where C::Point: Clone {}
