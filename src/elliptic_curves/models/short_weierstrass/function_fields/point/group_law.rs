use crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField;
use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    models::short_weierstrass::group_law_core::{
        ShortWeierstrassFormulaOps, ShortWeierstrassFormulaRunner,
    },
    short_weierstrass::function_fields::{
        ShortWeierstrassFunction, ShortWeierstrassFunctionFieldPoint,
    },
};
use crate::fields::rational_function_field::RationalFunction;
use crate::fields::traits::*;

pub(super) struct FunctionFieldOps<F: Field> {
    curve: ShortWeierstrassCurve<F>,
}

impl<F: Field> FunctionFieldOps<F> {
    fn new(curve: &ShortWeierstrassCurve<F>) -> Self {
        Self {
            curve: curve.clone(),
        }
    }
}

impl<F: Field> ShortWeierstrassFormulaOps for FunctionFieldOps<F> {
    type Coord = ShortWeierstrassFunction<F>;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        left.add(right)
    }

    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        left.sub(right)
    }

    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        left.mul(right)
    }

    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        value.inverse()
    }

    fn lift_i64(&self, value: i64) -> Self::Coord {
        ShortWeierstrassFunction::from_rational_function(
            self.curve.clone(),
            RationalFunction::constant(F::from_i64(value)),
        )
    }

    fn is_zero(&self, value: &Self::Coord) -> bool {
        value.is_zero()
    }

    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool {
        left == right
    }
}

impl<F: Field> ShortWeierstrassFunctionField<F> {
    #[cfg(test)]
    pub(crate) fn neg_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        point.ensure_point_lives_on(self.curve())?;
        match point {
            ShortWeierstrassFunctionFieldPoint::Infinity => {
                Ok(ShortWeierstrassFunctionFieldPoint::Infinity)
            }
            ShortWeierstrassFunctionFieldPoint::Affine { x, y } => {
                ShortWeierstrassFunctionFieldPoint::affine(x.clone(), y.neg())
            }
        }
    }

    pub(crate) fn add_points(
        &self,
        left: &ShortWeierstrassFunctionFieldPoint<F>,
        right: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        left.ensure_point_lives_on(self.curve())?;
        right.ensure_point_lives_on(self.curve())?;
        let ops = FunctionFieldOps::new(self.curve());
        let curve_a = ShortWeierstrassFunction::from_rational_function(
            self.curve().clone(),
            RationalFunction::constant(self.curve().a().clone()),
        );
        let runner = ShortWeierstrassFormulaRunner::new(&ops, &curve_a);
        let formula_point =
            runner.add_points(&left.to_formula_point(), &right.to_formula_point())?;
        ShortWeierstrassFunctionFieldPoint::from_formula_point(formula_point)
    }

    #[cfg(test)]
    pub(crate) fn double_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        point.ensure_point_lives_on(self.curve())?;
        let ops = FunctionFieldOps::new(self.curve());
        let curve_a = ShortWeierstrassFunction::from_rational_function(
            self.curve().clone(),
            RationalFunction::constant(self.curve().a().clone()),
        );
        let runner = ShortWeierstrassFormulaRunner::new(&ops, &curve_a);
        let formula_point = runner.double_point(&point.to_formula_point())?;
        ShortWeierstrassFunctionFieldPoint::from_formula_point(formula_point)
    }

    pub(crate) fn mul_scalar_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
        scalar: u64,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        point.ensure_point_lives_on(self.curve())?;
        let ops = FunctionFieldOps::new(self.curve());
        let curve_a = ShortWeierstrassFunction::from_rational_function(
            self.curve().clone(),
            RationalFunction::constant(self.curve().a().clone()),
        );
        let runner = ShortWeierstrassFormulaRunner::new(&ops, &curve_a);
        let formula_point = runner.mul_point(&point.to_formula_point(), scalar)?;
        ShortWeierstrassFunctionFieldPoint::from_formula_point(formula_point)
    }
}
