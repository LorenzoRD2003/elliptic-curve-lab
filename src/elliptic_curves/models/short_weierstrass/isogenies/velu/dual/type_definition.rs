use std::hash::Hash;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::{isogenies::VeluIsogeny, isomorphisms::ShortWeierstrassIsomorphism},
    traits::{CurveIsomorphism, CurveModel},
};
use crate::fields::traits::Field;
use crate::isogenies::{
    error::IsogenyError,
    kernel::KernelDescription,
    traits::{DegreeFactorizedIsogeny, Isogeny},
};

/// Exhaustively searched dual of a short-Weierstrass Vélu isogeny.
pub struct DualVeluIsogeny<F: Field> {
    velu_part: VeluIsogeny<ShortWeierstrassCurve<F>>,
    codomain_to_original: ShortWeierstrassIsomorphism<F>,
}

impl<F: Field> DualVeluIsogeny<F> {
    pub fn new(
        velu_part: VeluIsogeny<ShortWeierstrassCurve<F>>,
        codomain_to_original: ShortWeierstrassIsomorphism<F>,
    ) -> Self {
        Self {
            velu_part,
            codomain_to_original,
        }
    }

    pub fn velu_part(&self) -> &VeluIsogeny<ShortWeierstrassCurve<F>> {
        &self.velu_part
    }

    pub fn codomain_to_original(&self) -> &ShortWeierstrassIsomorphism<F> {
        &self.codomain_to_original
    }
}

impl<F: Field + Clone> Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for DualVeluIsogeny<F>
where
    F::Elem: Clone + Eq + Hash,
{
    fn domain(&self) -> &ShortWeierstrassCurve<F> {
        self.velu_part.domain()
    }

    fn codomain(&self) -> &ShortWeierstrassCurve<F> {
        self.codomain_to_original.codomain()
    }

    fn degree(&self) -> usize {
        self.velu_part.degree()
    }

    fn evaluate(
        &self,
        point: &<ShortWeierstrassCurve<F> as CurveModel>::Point,
    ) -> Result<<ShortWeierstrassCurve<F> as CurveModel>::Point, IsogenyError> {
        let mid = self.velu_part.evaluate(point)?;
        self.codomain_to_original.evaluate(&mid).map_err(Into::into)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        self.velu_part.kernel_description()
    }
}

impl<F: Field + Clone> DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for DualVeluIsogeny<F>
where
    F::Elem: Clone + Eq + Hash,
{
    fn separable_degree(&self) -> u128 {
        self.degree() as u128
    }

    fn inseparable_degree(&self) -> u128 {
        1
    }
}
