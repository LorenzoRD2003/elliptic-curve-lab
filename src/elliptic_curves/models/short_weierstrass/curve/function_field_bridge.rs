use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::rational_function_field::RationalFunction;
use crate::fields::traits::*;

impl<F: Field> ShortWeierstrassCurve<F> {
    pub(crate) fn function_field_curve_rhs_rational_function(&self) -> RationalFunction<F> {
        RationalFunction::from_polynomial(self.to_cubic())
    }

    pub(crate) fn evaluate_curve_rhs_at_function_x(
        &self,
        x_value: &ShortWeierstrassFunction<F>,
    ) -> Result<ShortWeierstrassFunction<F>, CurveError> {
        let x_squared = x_value.mul(x_value)?;
        let x_cubed = x_squared.mul(x_value)?;
        let a_times_x = x_value.mul(&ShortWeierstrassFunction::from_rational_function(
            x_value.curve().clone(),
            RationalFunction::<F>::constant(self.a().clone()),
        ))?;
        let constant_term = ShortWeierstrassFunction::from_rational_function(
            x_value.curve().clone(),
            RationalFunction::<F>::constant(self.b().clone()),
        );

        x_cubed.add(&a_times_x)?.add(&constant_term)
    }
}
