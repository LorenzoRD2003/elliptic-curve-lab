use core::marker::PhantomData;

use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve,
    models::short_weierstrass::group_law_core::{
        ShortWeierstrassFormulaOps, ShortWeierstrassFormulaPoint, ShortWeierstrassFormulaRunner,
    },
    traits::{CurveModel, GroupCurveModel},
};
use crate::fields::traits::Field;

struct BaseFieldOps<F: Field>(PhantomData<F>);

impl<F: Field> ShortWeierstrassFormulaOps for BaseFieldOps<F> {
    type Coord = F::Elem;

    fn add(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(F::add(left, right))
    }

    fn sub(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(F::sub(left, right))
    }

    fn mul(&self, left: &Self::Coord, right: &Self::Coord) -> Result<Self::Coord, CurveError> {
        Ok(F::mul(left, right))
    }

    fn inv(&self, value: &Self::Coord) -> Result<Self::Coord, CurveError> {
        F::inv(value).ok_or(CurveError::PointNotOnCurve)
    }

    fn lift_i64(&self, value: i64) -> Self::Coord {
        F::from_i64(value)
    }

    fn is_zero(&self, value: &Self::Coord) -> bool {
        F::is_zero(value)
    }

    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool {
        F::eq(left, right)
    }
}

impl<F: Field> ShortWeierstrassCurve<F> {
    fn formula_point_from_affine(
        &self,
        point: &AffinePoint<F>,
    ) -> ShortWeierstrassFormulaPoint<F::Elem> {
        match point {
            AffinePoint::Infinity => ShortWeierstrassFormulaPoint::Infinity,
            AffinePoint::Finite { x, y } => ShortWeierstrassFormulaPoint::Affine {
                x: x.clone(),
                y: y.clone(),
            },
        }
    }

    fn affine_point_from_formula(
        &self,
        point: ShortWeierstrassFormulaPoint<F::Elem>,
    ) -> AffinePoint<F> {
        match point {
            ShortWeierstrassFormulaPoint::Infinity => self.identity(),
            ShortWeierstrassFormulaPoint::Affine { x, y } => self.unchecked_point(x, y),
        }
    }

    /// Adds two curve points under the assumption that both are already valid.
    fn add_unchecked(
        &self,
        left: &AffinePoint<F>,
        right: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        let ops = BaseFieldOps::<F>(PhantomData);
        let runner = ShortWeierstrassFormulaRunner::new(&ops, self.a());
        let formula_point = runner.add_points(
            &self.formula_point_from_affine(left),
            &self.formula_point_from_affine(right),
        )?;
        let point = self.affine_point_from_formula(formula_point);
        debug_assert!(self.contains(&point));
        Ok(point)
    }

    /// Doubles a curve point under the assumption that it is already valid.
    fn double_unchecked(&self, point: &AffinePoint<F>) -> Result<AffinePoint<F>, CurveError> {
        let ops = BaseFieldOps::<F>(PhantomData);
        let runner = ShortWeierstrassFormulaRunner::new(&ops, self.a());
        let formula_point = runner.double_point(&self.formula_point_from_affine(point))?;
        let doubled = self.affine_point_from_formula(formula_point);
        debug_assert!(self.contains(&doubled));
        Ok(doubled)
    }

    /// Multiplies a valid curve point by a non-negative integer without
    /// repeating curve-membership checks along the way.
    fn mul_scalar_unchecked(
        &self,
        point: &AffinePoint<F>,
        scalar: u64,
    ) -> Result<AffinePoint<F>, CurveError> {
        let ops = BaseFieldOps::<F>(PhantomData);
        let runner = ShortWeierstrassFormulaRunner::new(&ops, self.a());
        let formula_point = runner.mul_point(&self.formula_point_from_affine(point), scalar)?;
        let result = self.affine_point_from_formula(formula_point);
        debug_assert!(self.contains(&result));
        Ok(result)
    }
}

impl<F: Field> GroupCurveModel for ShortWeierstrassCurve<F> {
    fn neg(&self, point: &Self::Point) -> Self::Point {
        point.neg()
    }

    fn add(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.add_unchecked(left, right)
    }

    fn sub(&self, left: &Self::Point, right: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(left) || !self.contains(right) {
            return Err(CurveError::PointNotOnCurve);
        }

        let negated = self.neg(right);
        self.add_unchecked(left, &negated)
    }

    fn double(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.double_unchecked(point)
    }

    fn mul_scalar(&self, point: &Self::Point, scalar: u64) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        self.mul_scalar_unchecked(point, scalar)
    }

    fn mul_scalar_signed(
        &self,
        point: &Self::Point,
        scalar: i64,
    ) -> Result<Self::Point, CurveError> {
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        if scalar < 0 {
            let negated = self.neg(point);
            self.mul_scalar_unchecked(&negated, scalar.unsigned_abs())
        } else {
            self.mul_scalar_unchecked(point, scalar as u64)
        }
    }
}
