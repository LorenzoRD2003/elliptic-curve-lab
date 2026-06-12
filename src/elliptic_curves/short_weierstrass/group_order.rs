use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::frobenius::{
    CharacterSumPointCount, FrobeniusTrace, GroupOrderReport, GroupOrderStrategy,
};
use crate::elliptic_curves::traits::FrobeniusTraceCurveModel;
use crate::fields::{
    EnumerableFiniteField, FiniteField, FiniteFieldDescriptor, QuadraticCharacterFiniteField,
    SqrtField,
};
use num_bigint::BigUint;

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    fn group_order_from_exponent_lower_bound(
        &self,
        exponent_lower_bound: BigUint,
        hasse_strategy: crate::elliptic_curves::HasseGroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        let verification = self.verify_exponent_lower_bound_against_group_order_report(
            exponent_lower_bound,
            self.group_order_by(hasse_strategy.as_group_order_strategy())?,
        );

        verification.verified_group_order().ok_or_else(|| {
            CurveError::UnverifiedGroupOrderFromExponentLowerBound {
                lower_bound: verification.exponent_lower_bound().clone(),
            }
        })?;

        Ok(GroupOrderReport::FromExponentLowerBound(Box::new(
            verification,
        )))
    }

    /// Computes `#E(F_q)` using one requested public strategy.
    ///
    /// This is the user-facing finite-field group-order entry point for the
    /// current short-Weierstrass model.
    ///
    /// Complexity:
    /// - [`GroupOrderStrategy::Exhaustive`]: dominated by full rational-point
    ///   enumeration
    /// - [`GroupOrderStrategy::QuadraticCharacter`]: `Θ(q)` right-hand-side
    ///   evaluations and quadratic-character queries
    /// - [`GroupOrderStrategy::Auto`]: currently the same as the
    ///   quadratic-character route
    pub fn group_order_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<GroupOrderReport, CurveError> {
        match strategy {
            GroupOrderStrategy::Auto | GroupOrderStrategy::QuadraticCharacter => self
                .group_order_by_quadratic_character()
                .map(GroupOrderReport::QuadraticCharacter),
            GroupOrderStrategy::Exhaustive => FrobeniusTraceCurveModel::frobenius_trace(self)
                .map(GroupOrderReport::ExhaustiveTrace),
            GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
                exponent_lower_bound,
                hasse_strategy,
            } => self.group_order_from_exponent_lower_bound(exponent_lower_bound, hasse_strategy),
        }
    }

    /// Internal short-Weierstrass-specific `Θ(q)` character-sum count.
    ///
    /// The public entry point for callers is [`Self::group_order_by`]. This
    /// helper stays crate-private so the public API has one primary counting
    /// door while still letting internal tests and visualizations exercise the
    /// specific route directly.
    ///
    /// Formula:
    /// `#E(F_q) = q + 1 + Σ_{x ∈ F_q} χ(x^3 + Ax + B)`.
    ///
    /// Complexity:
    /// `Θ(q)` evaluations of `x^3 + Ax + B` and `Θ(q)` quadratic-character
    /// queries over represented field elements.
    pub(crate) fn group_order_by_quadratic_character(
        &self,
    ) -> Result<CharacterSumPointCount, CurveError> {
        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: F::characteristic(),
                extension_degree: F::extension_degree().get(),
            })?;

        let mut character_sum = 0i128;
        for x in F::elements() {
            let rhs = self.rhs_value(&x);
            let value = F::quadratic_character_of(&rhs).map_err(|_| {
                CurveError::UnsupportedCharacterSumPointCount {
                    characteristic: F::characteristic(),
                    extension_degree: F::extension_degree().get(),
                }
            })?;
            character_sum += value.as_i128();
        }

        CharacterSumPointCount::new(base_field, character_sum)
    }

    /// Recovers the Frobenius trace through one requested group-order
    /// strategy.
    pub fn frobenius_trace_by(
        &self,
        strategy: GroupOrderStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.group_order_by(strategy)?.to_frobenius_trace()
    }
}
