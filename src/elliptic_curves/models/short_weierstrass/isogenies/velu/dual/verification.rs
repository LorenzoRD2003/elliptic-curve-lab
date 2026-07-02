use crate::fields::traits::*;
use std::hash::Hash;

use crate::elliptic_curves::short_weierstrass::isogenies::VeluIsogeny;
use crate::elliptic_curves::{ShortWeierstrassCurve, traits::EnumerableCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    error::{DualIsogenyError, IsogenyError},
    scalar_multiplication::ScalarMultiplicationIsogeny,
    traits::Isogeny,
};

use super::type_definition::DualVeluIsogeny;

impl<F: Field + EnumerableFiniteField + SqrtField + Clone> DualVeluIsogeny<F>
where
    F::Elem: Clone + Eq + Hash,
{
    /// Verifies exhaustively that `φ̂ ∘ φ = [n]_E` on rational points of the domain.
    pub fn verify_left_dual_relation(
        &self,
        phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    ) -> Result<(), IsogenyError> {
        let scalar = ScalarMultiplicationIsogeny::new(phi.domain().clone(), phi.degree())?;

        for point in phi.domain().points() {
            let left = self.evaluate(&phi.evaluate(&point)?)?;
            let right = scalar.evaluate(&point)?;
            if left != right {
                return Err(DualIsogenyError::DualRelationViolation.into());
            }
        }

        Ok(())
    }

    /// Verifies exhaustively that `φ ∘ φ̂ = [n]_{E'}` on rational points of the codomain.
    pub fn verify_right_dual_relation(
        &self,
        phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    ) -> Result<(), IsogenyError> {
        let scalar = ScalarMultiplicationIsogeny::new(phi.codomain().clone(), phi.degree())?;

        for point in phi.codomain().points() {
            let left = phi.evaluate(&self.evaluate(&point)?)?;
            let right = scalar.evaluate(&point)?;
            if left != right {
                return Err(DualIsogenyError::DualRelationViolation.into());
            }
        }

        Ok(())
    }

    /// Verifies both duality identities exhaustively on rational points.
    pub fn verify_dual_relations(
        &self,
        phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
    ) -> Result<(), IsogenyError> {
        self.verify_left_dual_relation(phi)?;
        self.verify_right_dual_relation(phi)
    }
}
