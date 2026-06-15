use crate::elliptic_curves::{AffinePoint, ShortWeierstrassCurve, frobenius::RelativeFrobenius};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    error::IsogenyError,
    kernel::{KernelDescription, NonReducedKernelDescription},
    traits::{DegreeFactorizedIsogeny, Isogeny},
};

/// Relative Frobenius isogeny.
///
/// `Frob_q: E -> E, (x, y) ↦ (x^q, y^q)` where `q = |F| = p^r`.
#[derive(Clone, Debug)]
pub struct RelativeFrobeniusIsogeny<F: FiniteField> {
    domain: ShortWeierstrassCurve<F>,
    codomain: ShortWeierstrassCurve<F>,
}

impl<F: FiniteField> RelativeFrobeniusIsogeny<F> {
    pub fn new(domain: ShortWeierstrassCurve<F>) -> Result<Self, IsogenyError> {
        let codomain = domain
            .frobenius_twist_power(F::extension_degree().get())
            .map_err(IsogenyError::Curve)?;

        Ok(Self { domain, codomain })
    }

    pub fn frobenius(&self) -> RelativeFrobenius {
        RelativeFrobenius::for_field::<F>(1)
    }
}

impl<F: FiniteField> DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for RelativeFrobeniusIsogeny<F>
{
    fn separable_degree(&self) -> u128 {
        1
    }

    fn inseparable_degree(&self) -> u128 {
        u128::from(F::characteristic())
            .checked_pow(F::extension_degree().get())
            .expect("relative Frobenius inseparable degree should fit in u128 in the educational setting")
    }
}

impl<F: FiniteField> Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    for RelativeFrobeniusIsogeny<F>
{
    fn domain(&self) -> &ShortWeierstrassCurve<F> {
        &self.domain
    }

    fn codomain(&self) -> &ShortWeierstrassCurve<F> {
        &self.codomain
    }

    fn degree(&self) -> usize {
        usize::try_from(self.total_degree())
            .expect("relative Frobenius degree should fit into usize in the educational setting")
    }

    fn evaluate(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, IsogenyError> {
        self.domain()
            .relative_frobenius_point(point)
            .map_err(IsogenyError::Curve)
    }

    fn kernel_description(&self) -> KernelDescription<ShortWeierstrassCurve<F>> {
        KernelDescription::NonReduced(NonReducedKernelDescription::new(
            self.degree(),
            "ker(Frob_q)",
        ))
    }
}
