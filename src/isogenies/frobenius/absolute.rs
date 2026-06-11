use crate::elliptic_curves::{
    AbsoluteFrobenius, CurveModel, ShortWeierstrassCurve, absolute_frobenius_power_point,
    frobenius_twist_power,
};
use crate::fields::FiniteField;
use crate::isogenies::{
    DegreeFactorizedIsogeny, Isogeny, IsogenyError, KernelDescription, NonReducedKernelDescription,
};

/// Absolute Frobenius isogeny
///
/// `Frob_p: E -> E^(p), (x, y) ↦ (x^p, y^p)`.
#[derive(Clone, Debug)]
pub struct AbsoluteFrobeniusIsogeny<F: FiniteField> {
    domain: ShortWeierstrassCurve<F>,
    codomain: ShortWeierstrassCurve<F>,
}

impl<F: FiniteField> AbsoluteFrobeniusIsogeny<F> {
    pub fn new(domain: ShortWeierstrassCurve<F>) -> Result<Self, IsogenyError> {
        let codomain = frobenius_twist_power(&domain, 1).map_err(IsogenyError::Curve)?;

        Ok(Self { domain, codomain })
    }

    pub fn frobenius(&self) -> AbsoluteFrobenius {
        AbsoluteFrobenius::for_field::<F>(1)
    }
}

impl<F: FiniteField> DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for AbsoluteFrobeniusIsogeny<F>
{
    fn separable_degree(&self) -> u128 {
        1
    }

    fn inseparable_degree(&self) -> u128 {
        u128::from(F::characteristic())
    }
}

impl<F: FiniteField> Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for AbsoluteFrobeniusIsogeny<F>
{
    fn domain(&self) -> &ShortWeierstrassCurve<F> {
        &self.domain
    }

    fn codomain(&self) -> &ShortWeierstrassCurve<F> {
        &self.codomain
    }

    fn degree(&self) -> usize {
        usize::try_from(self.total_degree())
            .expect("absolute Frobenius degree should fit into usize in the educational setting")
    }

    fn evaluate(
        &self,
        point: &<ShortWeierstrassCurve<F> as CurveModel>::Point,
    ) -> Result<<ShortWeierstrassCurve<F> as CurveModel>::Point, IsogenyError> {
        absolute_frobenius_power_point(self.domain(), point, 1).map_err(IsogenyError::Curve)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        KernelDescription::NonReduced(NonReducedKernelDescription::new(
            self.degree(),
            "ker(Frob_p)",
        ))
    }
}
