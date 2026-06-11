use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassFunction};
use crate::fields::FiniteField;
use crate::isogenies::frobenius::shared::{x_pullback_from_power, y_pullback_from_power};
use crate::isogenies::{
    DegreeFactorizedIsogeny, DifferentialPullbackReport, Isogeny, IsogenyError,
    ShortWeierstrassFunctionFieldMap,
};

/// Shared surface for short-Weierstrass Frobenius isogenies.
pub trait FrobeniusLikeIsogeny<F: FiniteField>:
    Isogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
    + DegreeFactorizedIsogeny<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>
{
    /// Returns the pullback of the codomain `x`-coordinate.
    fn x_pullback(&self) -> ShortWeierstrassFunction<F> {
        x_pullback_from_power(self.domain(), self.inseparable_degree())
    }

    /// Returns the pullback of the codomain `y`-coordinate.
    fn y_pullback(&self) -> ShortWeierstrassFunction<F> {
        y_pullback_from_power(self.domain(), self.inseparable_degree())
    }

    /// Returns the Frobenius pullback as a short-Weierstrass function-field map.
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

    /// Returns the differential pullback report.
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
