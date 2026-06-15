use crate::elliptic_curves::{
    CurveError, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{rational_function_field::RationalFunction, traits::Field};
use crate::polynomials::DensePolynomial;

impl<F: Field> ShortWeierstrassFunction<F> {
    pub(crate) fn evaluate_polynomial_at_function_x(
        polynomial: &DensePolynomial<F>,
        x_value: &Self,
    ) -> Result<Self, CurveError> {
        let mut accumulator = Self::zero(x_value.curve().clone());

        for coefficient in polynomial.coefficients().iter().rev() {
            accumulator = accumulator.mul(x_value)?;
            let coefficient_term = Self::from_rational_function(
                x_value.curve().clone(),
                RationalFunction::constant(coefficient.clone()),
            );
            accumulator = accumulator.add(&coefficient_term)?;
        }

        Ok(accumulator)
    }

    pub(crate) fn substitute_rational_function_at_function_x(
        function: &RationalFunction<F>,
        x_value: &Self,
    ) -> Result<Self, CurveError> {
        let numerator = Self::evaluate_polynomial_at_function_x(function.numerator(), x_value)?;
        let denominator = Self::evaluate_polynomial_at_function_x(function.denominator(), x_value)?;
        numerator.div(&denominator)
    }
}
