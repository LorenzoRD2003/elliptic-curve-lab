use crate::elliptic_curves::{
    CurveError, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{rational_function_field::RationalFunction, traits::Field};

impl<F: Field> ShortWeierstrassFunction<F> {
    pub fn conjugate(&self) -> Self {
        Self::new(
            self.curve().clone(),
            self.a_part().clone(),
            self.b_part().neg(),
        )
    }

    pub fn norm(&self) -> RationalFunction<F> {
        self.a_part().mul(self.a_part()).sub(
            &self
                .curve()
                .function_field_curve_rhs_rational_function()
                .mul(&self.b_part().mul(self.b_part())),
        )
    }

    pub fn add(&self, rhs: &Self) -> Result<Self, CurveError> {
        self.ensure_same_curve(rhs)?;
        Ok(Self::new(
            self.curve().clone(),
            self.a_part().add(rhs.a_part()),
            self.b_part().add(rhs.b_part()),
        ))
    }

    pub fn neg(&self) -> Self {
        Self::new(
            self.curve().clone(),
            self.a_part().neg(),
            self.b_part().neg(),
        )
    }

    pub fn sub(&self, rhs: &Self) -> Result<Self, CurveError> {
        self.ensure_same_curve(rhs)?;
        Ok(Self::new(
            self.curve().clone(),
            self.a_part().sub(rhs.a_part()),
            self.b_part().sub(rhs.b_part()),
        ))
    }

    pub fn mul(&self, rhs: &Self) -> Result<Self, CurveError> {
        self.ensure_same_curve(rhs)?;
        let (a, b, c, d) = (self.a_part(), self.b_part(), rhs.a_part(), rhs.b_part());

        let f = self.curve().function_field_curve_rhs_rational_function();
        let prod_a_part = a.mul(c).add(&f.mul(&b.mul(d)));
        let prod_b_part = a.mul(d).add(&b.mul(c));

        Ok(Self::new(self.curve().clone(), prod_a_part, prod_b_part))
    }

    pub fn inverse(&self) -> Result<Self, CurveError> {
        let norm = self.norm();
        if norm.is_zero() {
            return Err(CurveError::NonInvertibleFunctionFieldElement);
        }

        let a_part = self
            .a_part()
            .div(&norm)
            .map_err(|_| CurveError::NonInvertibleFunctionFieldElement)?;
        let b_part = self
            .b_part()
            .neg()
            .div(&norm)
            .map_err(|_| CurveError::NonInvertibleFunctionFieldElement)?;

        Ok(Self::new(self.curve().clone(), a_part, b_part))
    }

    pub fn div(&self, rhs: &Self) -> Result<Self, CurveError> {
        self.ensure_same_curve(rhs)?;
        self.mul(&rhs.inverse()?)
    }

    pub fn derivative(&self) -> Self {
        let f = self.curve().function_field_curve_rhs_rational_function();
        let two = RationalFunction::<F>::constant(F::from_i64(2));
        let correction = f
            .derivative()
            .mul(self.b_part())
            .div(&two.mul(&f))
            .expect("validated short-Weierstrass curve has non-zero cubic right-hand side in the rational function field");

        let a_part = self.a_part().derivative();
        let b_part = self.b_part().derivative().add(&correction);

        Self::new(self.curve().clone(), a_part, b_part)
    }
}
