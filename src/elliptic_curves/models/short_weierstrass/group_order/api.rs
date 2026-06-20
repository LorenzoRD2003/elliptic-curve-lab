use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        FrobeniusTrace,
        group_order::{
            FiniteFieldGroupOrderStrategy, GroupOrderReport, SchoofGroupOrderSummary,
            SmallFieldGroupOrderStrategy, SmallFieldSampledGroupOrderStrategy,
        },
    },
    short_weierstrass::division_polynomials::DivisionPolynomialError,
    traits::{FrobeniusTraceCurveModel, PointIndexSampler},
};
use crate::fields::traits::{
    EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField,
};

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    pub(crate) fn group_order_by_finite_field_strategy(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        match strategy {
            FiniteFieldGroupOrderStrategy::Auto | FiniteFieldGroupOrderStrategy::Schoof => self
                .schoof_group_order()
                .map_err(curve_error_from_schoof_division_polynomial_error)
                .and_then(|report| {
                    SchoofGroupOrderSummary::from_detailed(&report)
                        .map(Box::new)
                        .map(GroupOrderReport::Schoof)
                }),
        }
    }

    /// Computes `#E(F_q)` through one finite-field-capable route.
    ///
    /// The current finite-field surface is intentionally the general one: it
    /// does not assume exhaustive enumeration, quadratic-character support, or
    /// a square-root backend. Its `Auto` policy is therefore the automatic
    /// Schoof route.
    pub fn group_order_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        self.group_order_by_finite_field_strategy(strategy)
    }

    /// Recovers the Frobenius trace through one finite-field-capable
    /// group-order route.
    pub fn frobenius_trace_by(
        &self,
        strategy: FiniteFieldGroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }
}

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(crate) fn group_order_by_small_field_strategy(
        &self,
        strategy: SmallFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        match strategy {
            SmallFieldGroupOrderStrategy::Auto
            | SmallFieldGroupOrderStrategy::QuadraticCharacter => self
                .group_order_by_quadratic_character()
                .map(GroupOrderReport::QuadraticCharacter),
            SmallFieldGroupOrderStrategy::Exhaustive => {
                FrobeniusTraceCurveModel::frobenius_trace(self)
                    .map(GroupOrderReport::ExhaustiveTrace)
            }
            SmallFieldGroupOrderStrategy::Schoof => {
                self.group_order_by_finite_field_strategy(FiniteFieldGroupOrderStrategy::Schoof)
            }
        }
    }

    /// Computes `#E(F_q)` through one route that is specific to small
    /// enumerable finite fields.
    ///
    /// This surface exists for educational routes such as exhaustive point
    /// counting and quadratic-character counting. Its `Auto` policy remains the
    /// current small-field default: the quadratic-character route.
    pub fn group_order_by_small_field(
        &self,
        strategy: SmallFieldGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        self.group_order_by_small_field_strategy(strategy)
    }

    /// Recovers the Frobenius trace through one small-field deterministic
    /// group-order route.
    pub fn frobenius_trace_by_small_field(
        &self,
        strategy: SmallFieldGroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by_small_field(strategy)?
            .to_frobenius_trace()
    }

    /// Computes `#E(F_q)` through one small-field route, including
    /// sampler-driven routes such as Mestre's prime-field algorithm.
    ///
    /// The additional `F::Elem: Hash` bound is needed only because the current
    /// Mestre implementation reuses BSGS-style Hasse search machinery whose
    /// internal tables key on represented points and coordinates.
    pub fn group_order_by_small_field_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: SmallFieldSampledGroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<GroupOrderReport, CurveError>
    where
        F::Elem: Hash,
    {
        match strategy {
            SmallFieldSampledGroupOrderStrategy::Auto
            | SmallFieldSampledGroupOrderStrategy::QuadraticCharacter => self
                .group_order_by_small_field_strategy(
                    SmallFieldGroupOrderStrategy::QuadraticCharacter,
                ),
            SmallFieldSampledGroupOrderStrategy::Exhaustive => {
                self.group_order_by_small_field_strategy(SmallFieldGroupOrderStrategy::Exhaustive)
            }
            SmallFieldSampledGroupOrderStrategy::Schoof => {
                self.group_order_by_finite_field_strategy(FiniteFieldGroupOrderStrategy::Schoof)
            }
            SmallFieldSampledGroupOrderStrategy::MestreFp(config) => {
                self.group_order_by_mestre_fp(config, sampler)
            }
        }
    }

    /// Recovers the Frobenius trace through one small-field route, including
    /// sampler-driven routes.
    pub fn frobenius_trace_by_small_field_with_sampler<S: PointIndexSampler>(
        &self,
        strategy: SmallFieldSampledGroupOrderStrategy,
        sampler: &mut S,
    ) -> Result<FrobeniusTrace, CurveError>
    where
        F::Elem: Hash,
    {
        self.group_order_by_small_field_with_sampler(strategy, sampler)?
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
