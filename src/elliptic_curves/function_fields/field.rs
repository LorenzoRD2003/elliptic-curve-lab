use crate::elliptic_curves::short_weierstrass::group_law_core::{
    ShortWeierstrassFormulaCoordOps, ShortWeierstrassFormulaPoint, add_formula_points,
    double_formula_point, mul_formula_point,
};
use crate::elliptic_curves::traits::CurveModel;
use crate::elliptic_curves::{
    AffinePoint, CurveError, ShortWeierstrassCurve, ShortWeierstrassFunction,
    ShortWeierstrassFunctionFieldPoint,
};
use crate::fields::{AmbientField, Field, RationalFunction};

/// Runtime ambient family for the function field of one concrete
/// short-Weierstrass curve.
///
/// Unlike the base-field trait [`Field`], this ambient object
/// depends on a specific curve chosen at runtime.
pub struct ShortWeierstrassFunctionField<F: Field> {
    curve: ShortWeierstrassCurve<F>,
}

struct FunctionFieldOps<F: Field> {
    curve: ShortWeierstrassCurve<F>,
}

impl<F: Field> ShortWeierstrassFormulaCoordOps for FunctionFieldOps<F> {
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
            RationalFunction::<F>::constant(F::from_i64(value)),
        )
    }

    fn is_zero(&self, value: &Self::Coord) -> bool {
        value.is_zero()
    }

    fn eq(&self, left: &Self::Coord, right: &Self::Coord) -> bool {
        left == right
    }
}

fn function_field_point_to_formula_point<F: Field>(
    point: &ShortWeierstrassFunctionFieldPoint<F>,
) -> ShortWeierstrassFormulaPoint<ShortWeierstrassFunction<F>> {
    match point {
        ShortWeierstrassFunctionFieldPoint::Infinity => ShortWeierstrassFormulaPoint::Infinity,
        ShortWeierstrassFunctionFieldPoint::Affine { x, y } => {
            ShortWeierstrassFormulaPoint::Affine {
                x: x.clone(),
                y: y.clone(),
            }
        }
    }
}

fn formula_point_to_function_field_point<F: Field>(
    point: ShortWeierstrassFormulaPoint<ShortWeierstrassFunction<F>>,
) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
    match point {
        ShortWeierstrassFormulaPoint::Infinity => Ok(ShortWeierstrassFunctionFieldPoint::Infinity),
        ShortWeierstrassFormulaPoint::Affine { x, y } => {
            ShortWeierstrassFunctionFieldPoint::affine(x, y)
        }
    }
}

impl<F: Field> ShortWeierstrassFunctionField<F> {
    /// Builds the ambient function field of a short Weierstrass curve.
    pub fn new(curve: ShortWeierstrassCurve<F>) -> Self {
        Self { curve }
    }

    /// Returns the ambient short-Weierstrass curve.
    pub fn curve(&self) -> &ShortWeierstrassCurve<F> {
        &self.curve
    }

    /// Returns the zero function.
    pub fn zero(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::zero(self.curve.clone())
    }

    /// Returns the constant function `1`.
    pub fn one(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::one(self.curve.clone())
    }

    /// Returns the coordinate function `x`.
    pub fn x(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::from_rational_function(
            self.curve.clone(),
            RationalFunction::<F>::indeterminate(),
        )
    }

    /// Returns the coordinate function `y`.
    pub fn y(&self) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::y(self.curve.clone())
    }

    /// Returns the generic affine point `(x, y)` of the current curve.
    ///
    /// This is the function-field incarnation of the tautological point on
    /// `E` whose coordinates are the distinguished generators of `F(E)`.
    ///
    /// Concretely, if `E: y^2 = x^3 + ax + b`, then this returns the affine
    /// function-field point `P_gen = (x, y) in E(F(E))`.
    ///
    /// This is the natural starting point for deriving pullback formulas of
    /// endomorphisms and isogenies: to compute `φ*(x')` and `φ*(y')`,
    /// one evaluates the target coordinate formulas at this point.
    pub fn generic_point(&self) -> ShortWeierstrassFunctionFieldPoint<F> {
        ShortWeierstrassFunctionFieldPoint::Affine {
            x: self.x(),
            y: self.y(),
        }
    }

    /// Embeds a base-field point of `E(F)` as a constant point of `E(F(E))`.
    pub(crate) fn embed_affine_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        if !self.curve.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(ShortWeierstrassFunctionFieldPoint::Infinity),
            AffinePoint::Finite { .. } => ShortWeierstrassFunctionFieldPoint::affine(
                ShortWeierstrassFunction::<F>::from_finite_point_coordinate(
                    self.curve.clone(),
                    point,
                    true,
                ),
                ShortWeierstrassFunction::<F>::from_finite_point_coordinate(
                    self.curve.clone(),
                    point,
                    false,
                ),
            ),
        }
    }

    /// Negates a function-field point.
    #[cfg(test)]
    pub(crate) fn neg_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.ensure_point_lives_on_this_curve(point)?;

        match point {
            ShortWeierstrassFunctionFieldPoint::Infinity => {
                Ok(ShortWeierstrassFunctionFieldPoint::Infinity)
            }
            ShortWeierstrassFunctionFieldPoint::Affine { x, y } => {
                ShortWeierstrassFunctionFieldPoint::affine(x.clone(), y.neg())
            }
        }
    }

    /// Adds two points of `E(F(E))` using the affine short-Weierstrass group law.
    pub(crate) fn add_points(
        &self,
        left: &ShortWeierstrassFunctionFieldPoint<F>,
        right: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.ensure_point_lives_on_this_curve(left)?;
        self.ensure_point_lives_on_this_curve(right)?;
        let ops = FunctionFieldOps {
            curve: self.curve.clone(),
        };
        let formula_point = add_formula_points(
            &ops,
            &ShortWeierstrassFunction::from_rational_function(
                self.curve.clone(),
                RationalFunction::<F>::constant(self.curve.a().clone()),
            ),
            &function_field_point_to_formula_point(left),
            &function_field_point_to_formula_point(right),
        )?;
        formula_point_to_function_field_point(formula_point)
    }

    /// Doubles one point of `E(F(E))` using the tangent formula.
    pub(crate) fn double_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.ensure_point_lives_on_this_curve(point)?;
        let ops = FunctionFieldOps {
            curve: self.curve.clone(),
        };
        let formula_point = double_formula_point(
            &ops,
            &ShortWeierstrassFunction::from_rational_function(
                self.curve.clone(),
                RationalFunction::<F>::constant(self.curve.a().clone()),
            ),
            &function_field_point_to_formula_point(point),
        )?;
        formula_point_to_function_field_point(formula_point)
    }

    /// Multiplies a function-field point by a non-negative integer using double-and-add.
    pub(crate) fn mul_scalar_point(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
        scalar: u64,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.ensure_point_lives_on_this_curve(point)?;
        let ops = FunctionFieldOps {
            curve: self.curve.clone(),
        };
        let formula_point = mul_formula_point(
            &ops,
            &ShortWeierstrassFunction::from_rational_function(
                self.curve.clone(),
                RationalFunction::<F>::constant(self.curve.a().clone()),
            ),
            &function_field_point_to_formula_point(point),
            scalar,
        )?;
        formula_point_to_function_field_point(formula_point)
    }

    /// Translates the generic point by one concrete affine point of the curve.
    ///
    /// If `Q` is a base-field point of `E`, this returns the function-field
    /// point `P_gen + Q`.
    ///
    /// Internally, this embeds `Q` as a constant point of `E(F(E))` and then
    /// applies [`Self::add_points`] to the generic point `P_gen = (x, y)`.
    ///
    /// Complexity: `Θ(1)` function-field additions, multiplications, and one
    /// inversion, since it performs one affine group-law addition.
    pub fn translate_generic_point_by_finite_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        let generic_point = self.generic_point();
        let constant_point = self.embed_affine_point(point)?;
        self.add_points(&generic_point, &constant_point)
    }

    /// Doubles the generic affine point inside `E(F(E))`.
    ///
    /// This returns the function-field point
    ///
    /// `[2]P_gen = (x([2]P), y([2]P))`
    ///
    /// using the usual tangent formulas on the short-Weierstrass model:
    ///
    /// - `λ = (3x^2 + a) / (2y)`,
    /// - `x([2]P) = λ^2 - 2x`,
    /// - `y([2]P) = λ(x - x([2]P)) - y`.
    ///
    /// Internally, this is just [`Self::double_point`] applied to the generic
    /// point.
    ///
    /// Complexity: `Θ(1)` function-field additions, multiplications, and one
    /// inversion, since it performs one affine tangent-formula doubling.
    pub fn double_generic_point(
        &self,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.double_point(&self.generic_point())
    }

    /// Returns `[n]P_gen` for the generic point `P_gen = (x, y)`.
    ///
    /// Complexity: `Θ(log n)` function-field additions/doublings under the
    /// current double-and-add backend.
    pub fn generic_point_multiple(
        &self,
        scalar: u64,
    ) -> Result<ShortWeierstrassFunctionFieldPoint<F>, CurveError> {
        self.mul_scalar_point(&self.generic_point(), scalar)
    }

    /// Embeds a rational function of `x` as `A(x) + y * 0`.
    pub fn from_rational_function(
        &self,
        function: RationalFunction<F>,
    ) -> ShortWeierstrassFunction<F> {
        ShortWeierstrassFunction::<F>::from_rational_function(self.curve.clone(), function)
    }

    fn ensure_point_lives_on_this_curve(
        &self,
        point: &ShortWeierstrassFunctionFieldPoint<F>,
    ) -> Result<(), CurveError> {
        match point {
            ShortWeierstrassFunctionFieldPoint::Infinity => Ok(()),
            ShortWeierstrassFunctionFieldPoint::Affine { x, y } => {
                if !F::eq(x.curve().a(), self.curve.a()) || !F::eq(x.curve().b(), self.curve.b()) {
                    return Err(CurveError::IncompatibleFunctionFieldCurves);
                }
                if !F::eq(y.curve().a(), self.curve.a()) || !F::eq(y.curve().b(), self.curve.b()) {
                    return Err(CurveError::IncompatibleFunctionFieldCurves);
                }

                let lhs = y.mul(y)?;
                let rhs = ShortWeierstrassFunction::evaluate_curve_rhs_in_x(&self.curve, x)?;
                if lhs == rhs {
                    Ok(())
                } else {
                    Err(CurveError::PointNotOnCurve)
                }
            }
        }
    }
}

impl<F: Field> Clone for ShortWeierstrassFunctionField<F> {
    fn clone(&self) -> Self {
        Self {
            curve: self.curve.clone(),
        }
    }
}

impl<F: Field> AmbientField for ShortWeierstrassFunctionField<F> {
    type Elem = ShortWeierstrassFunction<F>;
    type Error = CurveError;

    fn zero(&self) -> Self::Elem {
        Self::zero(self)
    }

    fn one(&self) -> Self::Elem {
        Self::one(self)
    }

    fn eq(&self, left: &Self::Elem, right: &Self::Elem) -> bool {
        left == right
    }

    fn add(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        left.add(right)
    }

    fn neg(&self, value: &Self::Elem) -> Self::Elem {
        value.neg()
    }

    fn mul(&self, left: &Self::Elem, right: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        left.mul(right)
    }

    fn inverse(&self, value: &Self::Elem) -> Result<Self::Elem, Self::Error> {
        value.inverse()
    }
}
