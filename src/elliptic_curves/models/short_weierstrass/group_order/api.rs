use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        group_order::{GroupOrderReport, GroupOrderStrategy, SchoofGroupOrderSummary},
    },
    short_weierstrass::division_polynomials::DivisionPolynomialError,
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
            GroupOrderStrategy::Schoof => self
                .schoof_group_order()
                .map_err(curve_error_from_schoof_division_polynomial_error)
                .and_then(|report| {
                    SchoofGroupOrderSummary::from_detailed(&report)
                        .map(Box::new)
                        .map(GroupOrderReport::Schoof)
                }),
            GroupOrderStrategy::MestreFp(_) => Err(CurveError::GroupOrderStrategyRequiresSampler {
                strategy: "MestreFp",
            }),
        }
    }

    /// Computes `#E(F_q)` using one requested public strategy.
    ///
    /// Note: this shared dispatcher currently lives on the same small-finite
    /// impl block as the exhaustive and quadratic-character routes. So even
    /// though the detailed Schoof implementation itself only needs
    /// `F: FiniteField`, the integrated `GroupOrderStrategy::Schoof` variant is
    /// still exposed here together with the stronger educational backend
    /// bounds needed by the other deterministic routes.
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

fn curve_error_from_schoof_division_polynomial_error(error: DivisionPolynomialError) -> CurveError {
    match error {
        DivisionPolynomialError::Curve(curve_error) => curve_error,
        DivisionPolynomialError::ZeroIndex
        | DivisionPolynomialError::UnsupportedIndex { .. }
        | DivisionPolynomialError::EvenIndexRequiresYFactor { .. }
        | DivisionPolynomialError::PointAtInfinityNotSupported => {
            CurveError::SchoofUnexpectedDivisionPolynomialModel
        }
        DivisionPolynomialError::Polynomial(error) => CurveError::SchoofPolynomialFailure { error },
        DivisionPolynomialError::FieldNotEnumerable => {
            CurveError::SchoofUnexpectedFieldEnumerationRequirement
        }
        DivisionPolynomialError::FieldHasNoSquareRootBackend => {
            CurveError::SchoofUnexpectedSquareRootRequirement
        }
    }
}
