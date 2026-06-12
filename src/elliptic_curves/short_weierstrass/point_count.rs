use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::frobenius::{
    CharacterSumPointCount, FrobeniusTrace, PointCountReport, PointCountStrategy,
};
use crate::elliptic_curves::traits::FrobeniusTraceCurveModel;
use crate::fields::{
    EnumerableFiniteField, FiniteField, FiniteFieldDescriptor, QuadraticCharacterFiniteField,
    SqrtField,
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    /// Counts `#E(F_q)` using one requested public strategy.
    ///
    /// This is the user-facing finite-field point-count entry point for the
    /// current short-Weierstrass model.
    ///
    /// Complexity:
    /// - [`PointCountStrategy::Exhaustive`]: dominated by full rational-point
    ///   enumeration
    /// - [`PointCountStrategy::QuadraticCharacter`]: `Θ(q)` right-hand-side
    ///   evaluations and quadratic-character queries
    /// - [`PointCountStrategy::Auto`]: currently the same as the
    ///   quadratic-character route
    pub fn count_points(&self, method: PointCountStrategy) -> Result<PointCountReport, CurveError> {
        match method {
            PointCountStrategy::Auto | PointCountStrategy::QuadraticCharacter => self
                .count_points_by_quadratic_character()
                .map(PointCountReport::QuadraticCharacter),
            PointCountStrategy::Exhaustive => FrobeniusTraceCurveModel::frobenius_trace(self)
                .map(PointCountReport::ExhaustiveTrace),
        }
    }

    /// Internal short-Weierstrass-specific `Θ(q)` character-sum count.
    ///
    /// The public entry point for callers is [`Self::count_points`]. This
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
    pub(crate) fn count_points_by_quadratic_character(
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

    /// Recovers the Frobenius trace through one requested counting strategy.
    pub fn frobenius_trace_by(
        &self,
        strategy: PointCountStrategy,
    ) -> Result<FrobeniusTrace, CurveError> {
        self.count_points(strategy)?.to_frobenius_trace()
    }
}
