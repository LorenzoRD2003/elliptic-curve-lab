use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        group_order::{GroupOrderReport, GroupOrderStrategy},
    },
    traits::{FrobeniusTraceCurveModel, PointIndexSampler},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(crate) fn group_order_by_deterministic_strategy(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        match strategy {
            GroupOrderStrategy::Auto | GroupOrderStrategy::QuadraticCharacter => self
                .group_order_by_quadratic_character()
                .map(GroupOrderReport::QuadraticCharacter),
            GroupOrderStrategy::Exhaustive => FrobeniusTraceCurveModel::frobenius_trace(self)
                .map(GroupOrderReport::ExhaustiveTrace),
            GroupOrderStrategy::MestreFp(_) => Err(CurveError::GroupOrderStrategyRequiresSampler {
                strategy: "MestreFp",
            }),
        }
    }

    /// Computes `#E(F_q)` using one requested public strategy.
    pub fn group_order_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        self.group_order_by_deterministic_strategy(strategy)
    }

    /// Computes `#E(F_q)` using one requested public strategy, including
    /// sampler-driven routes such as Mestre's prime-field algorithm.
    ///
    /// The additional `F::Elem: Hash` bound is needed only because the current
    /// Mestre implementation reuses BSGS-style Hasse search machinery whose
    /// internal tables key on represented points and coordinates.
    pub fn group_order_by_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: GroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<GroupOrderReport, CurveError>
    where
        F::Elem: Hash,
    {
        match strategy {
            GroupOrderStrategy::MestreFp(config) => self.group_order_by_mestre_fp(config, sampler),
            other => self.group_order_by_deterministic_strategy(other),
        }
    }

    /// Recovers the Frobenius trace through one requested group-order
    /// strategy for deterministic routes.
    pub fn frobenius_trace_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }

    /// Recovers the Frobenius trace through one requested group-order
    /// strategy, including sampler-driven routes.
    ///
    /// Like [`Self::group_order_by_with_sampler`], this sampler-aware variant
    /// exists so algorithms such as Mestre can consume externally chosen point
    /// samples without inventing hidden randomness defaults.
    pub fn frobenius_trace_by_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: GroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<FrobeniusTrace, CurveError>
    where
        F::Elem: Hash,
    {
        self.group_order_by_with_sampler(strategy, sampler)?
            .to_frobenius_trace()
    }
}
