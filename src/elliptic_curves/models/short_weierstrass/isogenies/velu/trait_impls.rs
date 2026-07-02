use crate::fields::traits::*;
use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve, short_weierstrass::isogenies::VeluIsogeny,
    traits::CurveModel,
};
use crate::isogenies::{
    error::IsogenyError,
    kernel::{KernelDescription, ReducedKernelDescription},
    traits::{DegreeFactorizedIsogeny, Isogeny},
};

impl<F: Field + Clone> Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F::Elem: Clone + Eq + Hash,
{
    fn domain(&self) -> &ShortWeierstrassCurve<F> {
        &self.domain
    }

    fn codomain(&self) -> &ShortWeierstrassCurve<F> {
        &self.codomain
    }

    fn degree(&self) -> usize {
        self.degree
    }

    fn evaluate(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, IsogenyError> {
        if !self.domain.contains(point) {
            return Err(CurveError::PointNotOnCurve.into());
        }

        if self.kernel.contains(point) {
            return Ok(self.codomain.identity());
        }

        let image = self.evaluate_non_kernel_point(point)?;
        debug_assert!(self.codomain.contains(&image));
        Ok(image)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        KernelDescription::Reduced(ReducedKernelDescription::RationalPointSubgroup(
            self.kernel.clone(),
        ))
    }
}

impl<F: Field + Clone> DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F::Elem: Clone + Eq + Hash,
{
    fn separable_degree(&self) -> num_bigint::BigUint {
        num_bigint::BigUint::from(self.degree)
    }

    fn inseparable_degree(&self) -> num_bigint::BigUint {
        num_bigint::BigUint::from(1u8)
    }
}
