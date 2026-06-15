use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, frobenius::character_sum::CharacterSumPointCount,
};
use crate::fields::{
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField},
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
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
}
