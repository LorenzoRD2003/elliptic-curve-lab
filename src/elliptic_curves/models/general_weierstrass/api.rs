use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        group_order::{
            FiniteFieldGroupOrderStrategy, GroupOrderReport, SmallFieldGroupOrderStrategy,
        },
    },
    group_algorithms::{shared_group_exponent_by, shared_point_order_by},
    short_weierstrass::{
        group_exponent::{GroupExponentReport, GroupExponentStrategy},
        point_order::{PointOrderReport, PointOrderStrategy},
    },
    traits::{
        CurveModelConversion, FiniteGroupCurveModel, FrobeniusTraceCurveModel, PointIndexSampler,
    },
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};
impl<F: FiniteField> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    /// Computes `#E(F_q)` through one finite-field-capable route.
    ///
    /// Current status:
    /// - `Auto` and `Schoof` are delegated to the short-Weierstrass companion
    ///   when the classical reduction exists
    /// - this wrapper is therefore unavailable in characteristics `2` and `3`
    ///   until the general model gains its own non-enumerative finite-field
    ///   counting route
    pub fn group_order_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        let conversion = self.conversion_to_short_weierstrass()?;
        conversion.target().group_order_by(strategy)
    }

    /// Recovers the Frobenius trace through one finite-field-capable
    /// group-order route.
    ///
    /// Current status:
    /// - delegated to the short-Weierstrass companion
    /// - unavailable in characteristics `2` and `3` for now
    pub fn frobenius_trace_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }
}

impl<F> GeneralWeierstrassCurve<F>
where
    F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField,
    F::Elem: Clone,
{
    /// Computes `#E(F_q)` through one route that is specific to small
    /// enumerable finite fields.
    ///
    /// Current status:
    /// - `Exhaustive` is native to `GeneralWeierstrassCurve<F>` through direct
    ///   point enumeration
    /// - `Auto`, `QuadraticCharacter`, and `Schoof` reuse the short companion
    ///   when available
    /// - `Auto` falls back to the native exhaustive route in characteristics
    ///   `2` and `3`, where the short reduction is unavailable
    pub fn group_order_by_small_field(
        &self,
        strategy: SmallFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError>
    where
        Self: FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
        match strategy {
            SmallFieldGroupOrderStrategy::Exhaustive => {
                FrobeniusTraceCurveModel::frobenius_trace(self)
                    .map(GroupOrderReport::ExhaustiveTrace)
            }
            SmallFieldGroupOrderStrategy::Auto
                if F::has_characteristic(2) || F::has_characteristic(3) =>
            {
                FrobeniusTraceCurveModel::frobenius_trace(self)
                    .map(GroupOrderReport::ExhaustiveTrace)
            }
            SmallFieldGroupOrderStrategy::Auto
            | SmallFieldGroupOrderStrategy::QuadraticCharacter
            | SmallFieldGroupOrderStrategy::Schoof => {
                let conversion = self.conversion_to_short_weierstrass()?;
                conversion.target().group_order_by_small_field(strategy)
            }
        }
    }

    /// Recovers the exact order of one point by one requested strategy.
    ///
    /// Current status:
    /// - `Exhaustive` is native
    /// - `FromKnownMultiple` is native through the shared big-scalar and
    ///   cyclic-primary helpers
    /// - `HasseIntervalNaive` is native once the chosen small-field
    ///   group-order route produces a Hasse interval
    pub fn point_order_by(
        &self,
        point: &AffinePoint<F>,
        strategy: PointOrderStrategy,
    ) -> Result<PointOrderReport<AffinePoint<F>>, CurveError>
    where
        Self: FiniteGroupCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
            + FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
        shared_point_order_by(self, point, strategy, |group_order_strategy| {
            self.group_order_by_small_field(group_order_strategy)
        })
    }

    /// Recovers or estimates `λ(E(F_q))` by one requested strategy.
    ///
    /// Current status:
    /// - `Exhaustive` is native
    /// - `RandomPoints` is native once the chosen point-order route is
    ///   available
    pub fn group_exponent_by<S: PointIndexSampler>(
        &self,
        strategy: GroupExponentStrategy,
        sampler: &mut S,
    ) -> Result<GroupExponentReport<AffinePoint<F>>, CurveError>
    where
        Self: FiniteGroupCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
            + FrobeniusTraceCurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>,
    {
        shared_group_exponent_by(self, strategy, sampler, |point, point_order_strategy| {
            self.point_order_by(point, point_order_strategy)
        })
    }
}
