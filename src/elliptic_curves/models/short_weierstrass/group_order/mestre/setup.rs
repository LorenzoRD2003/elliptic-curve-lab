use std::hash::Hash;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::{ShortWeierstrassQuadraticTwist, TwistKind},
};
use crate::fields::{
    finite_field_descriptor::FiniteFieldDescriptor,
    traits::{EnumerableFiniteField, FiniteField, QuadraticCharacterFiniteField, SqrtField},
};

impl<F: EnumerableFiniteField + FiniteField + QuadraticCharacterFiniteField + SqrtField>
    ShortWeierstrassCurve<F>
{
    pub(super) fn validate_mestre_prime_field(
        &self,
    ) -> Result<(FiniteFieldDescriptor, u128), CurveError> {
        if F::extension_degree().get() != 1 {
            return Err(CurveError::MestreRequiresPrimeField {
                extension_degree: F::extension_degree().get(),
            });
        }

        if F::characteristic() <= 229 {
            return Err(CurveError::MestrePrimeTooSmall {
                characteristic: F::characteristic(),
            });
        }

        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: F::characteristic(),
                extension_degree: F::extension_degree().get(),
            })?;
        let prime =
            base_field
                .cardinality()
                .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                    characteristic: F::characteristic(),
                    extension_degree: F::extension_degree().get(),
                })?;

        Ok((base_field, prime))
    }

    pub(super) fn select_genuine_quadratic_twist_for_mestre(
        &self,
    ) -> Result<ShortWeierstrassCurve<F>, CurveError>
    where
        F::Elem: Hash,
    {
        for candidate in F::elements() {
            if F::is_zero(&candidate) {
                continue;
            }

            let Ok(package) = ShortWeierstrassQuadraticTwist::new(self.clone(), candidate) else {
                continue;
            };
            if package.kind() == TwistKind::Quadratic {
                return Ok(package.twist().clone());
            }
        }

        Err(CurveError::MestreQuadraticTwistUnavailable)
    }
}
