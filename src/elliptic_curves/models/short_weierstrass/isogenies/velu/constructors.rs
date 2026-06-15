use std::collections::HashSet;
use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    short_weierstrass::isogenies::{VeluIsogeny, velu::VeluKernelData},
    traits::CurveModel,
};
use crate::fields::traits::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::{error::IsogenyError, kernel::IsogenyKernel};

impl<F: Field> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F::Elem: Clone + Eq + Hash,
{
    pub fn from_points(
        domain: ShortWeierstrassCurve<F>,
        kernel_points: HashSet<AffinePoint<F>>,
    ) -> Result<Self, IsogenyError> {
        let kernel = IsogenyKernel::new(&domain, kernel_points)?;
        Self::from_kernel(domain, kernel)
    }

    pub fn from_generator(
        domain: ShortWeierstrassCurve<F>,
        generator: AffinePoint<F>,
    ) -> Result<Self, IsogenyError>
    where
        F: EnumerableFiniteField + SqrtField,
        F::Elem: PartialEq,
    {
        let kernel = IsogenyKernel::cyclic(&domain, &generator)?;
        Self::from_kernel(domain, kernel)
    }

    pub(crate) fn from_kernel(
        domain: ShortWeierstrassCurve<F>,
        kernel: IsogenyKernel<ShortWeierstrassCurve<F>>,
    ) -> Result<Self, IsogenyError> {
        let codomain = VeluKernelData::from_kernel(&domain, &kernel).codomain_curve(&domain)?;
        Ok(Self::from_parts(domain, codomain, kernel))
    }

    pub(crate) fn require_non_kernel_finite_point<'a>(
        &self,
        point: &'a AffinePoint<F>,
    ) -> Result<Option<&'a AffinePoint<F>>, IsogenyError> {
        if !self.domain.contains(point) {
            return Err(CurveError::PointNotOnCurve.into());
        }

        if self.kernel.contains(point) {
            return Ok(None);
        }

        Ok((AffinePoint::finite_coordinates(point).is_some()).then_some(point))
    }
}
