use std::hash::Hash;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isogenies::{VeluIsogeny, velu::dual::type_definition::DualVeluIsogeny},
    traits::FiniteGroupCurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::{
    error::{DualIsogenyError, IsogenyError},
    traits::Isogeny,
};

impl<F> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash,
{
    fn candidate_satisfies_dual_relations_exhaustively(
        &self,
        candidate_dual: &DualVeluIsogeny<F>,
    ) -> Result<bool, IsogenyError> {
        match candidate_dual.verify_dual_relations(self) {
            Ok(()) => Ok(true),
            Err(IsogenyError::Dual(DualIsogenyError::DualRelationViolation)) => Ok(false),
            Err(other) => Err(other),
        }
    }

    /// Searches exhaustively for a dual isogeny on a small finite curve.
    pub fn find_dual_exhaustively(&self) -> Result<DualVeluIsogeny<F>, IsogenyError> {
        let degree = self.degree();
        let mut saw_isomorphic_candidate = false;

        for generator in self.codomain().points_of_order(degree) {
            let psi = VeluIsogeny::from_generator(self.codomain().clone(), generator)?;
            if psi.degree() != degree {
                return Err(DualIsogenyError::DegreeMismatch.into());
            }

            let alphas = psi.codomain().exhaustive_isomorphisms_to(self.domain());
            if alphas.is_empty() {
                continue;
            }
            saw_isomorphic_candidate = true;

            for alpha in alphas {
                let candidate_dual = DualVeluIsogeny::new(psi.clone(), alpha);

                if self.candidate_satisfies_dual_relations_exhaustively(&candidate_dual)? {
                    if candidate_dual.degree() != degree {
                        return Err(DualIsogenyError::DegreeMismatch.into());
                    }
                    return Ok(candidate_dual);
                }
            }
        }

        if saw_isomorphic_candidate {
            Err(DualIsogenyError::DualRelationViolation.into())
        } else {
            Err(DualIsogenyError::DualNotFound.into())
        }
    }
}
