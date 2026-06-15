use crate::elliptic_curves::short_weierstrass::{
    ShortWeierstrassCurve,
    function_fields::ShortWeierstrassFunction,
    isogenies::frobenius::shared::{x_pullback_from_power, y_pullback_from_power},
    isogenies::function_field_maps::{
        DifferentialPullbackReport, ShortWeierstrassFunctionFieldMap,
    },
};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    error::IsogenyError,
    traits::{DegreeFactorizedIsogeny, Isogeny},
};

/// Shared surface for short-Weierstrass Frobenius isogenies.
pub trait FrobeniusLikeIsogeny<F: FiniteField>:
    Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    + DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
{
    fn x_pullback(&self) -> ShortWeierstrassFunction<F> {
        x_pullback_from_power(self.domain(), self.inseparable_degree())
    }

    fn y_pullback(&self) -> ShortWeierstrassFunction<F> {
        y_pullback_from_power(self.domain(), self.inseparable_degree())
    }

    fn as_function_field_map(&self) -> ShortWeierstrassFunctionFieldMap<F>
    where
        F::Elem: PartialEq,
    {
        ShortWeierstrassFunctionFieldMap::new(
            self.domain().clone(),
            self.codomain().clone(),
            self.x_pullback(),
            self.y_pullback(),
        )
        .expect("Frobenius pullbacks should satisfy the codomain equation")
    }

    fn differential_pullback_report(&self) -> Result<DifferentialPullbackReport<F>, IsogenyError>
    where
        F::Elem: PartialEq,
    {
        self.as_function_field_map().differential_pullback_report()
    }
}

impl<F, T> FrobeniusLikeIsogeny<F> for T
where
    F: FiniteField,
    T: Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
        + DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>,
{
}
