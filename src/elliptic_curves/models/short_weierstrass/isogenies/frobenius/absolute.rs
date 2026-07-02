use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve, frobenius::AbsoluteFrobenius};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    error::IsogenyError,
    kernel::{KernelDescription, NonReducedKernelDescription},
    traits::{DegreeFactorizedIsogeny, Isogeny},
};
use num_bigint::BigUint;
use num_traits::ToPrimitive;

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
        let codomain = domain
            .frobenius_twist_power(1)
            .map_err(IsogenyError::Curve)?;

        Ok(Self { domain, codomain })
    }

    pub fn frobenius(&self) -> AbsoluteFrobenius {
        AbsoluteFrobenius::for_field::<F>(1)
    }
}

impl<F: FiniteField> DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for AbsoluteFrobeniusIsogeny<F>
{
    fn separable_degree(&self) -> BigUint {
        BigUint::from(1u8)
    }

    fn inseparable_degree(&self) -> BigUint {
        F::characteristic()
            .to_positive_biguint()
            .expect("finite fields have positive characteristic")
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
        self.total_degree()
            .to_usize()
            .expect("absolute Frobenius degree should fit into usize in the educational setting")
    }

    fn evaluate(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, IsogenyError> {
        self.domain()
            .absolute_frobenius_power_point(point, 1)
            .map_err(IsogenyError::Curve)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        KernelDescription::NonReduced(NonReducedKernelDescription::new(
            self.degree(),
            "ker(Frob_p)",
        ))
    }
}
