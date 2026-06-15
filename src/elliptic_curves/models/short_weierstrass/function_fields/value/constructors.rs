use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::{rational_function_field::RationalFunction, traits::Field};

impl<F: Field> ShortWeierstrassFunction<F> {
    /// Builds the raw function-field element `A(x) + y B(x)` on a fixed
    /// short-Weierstrass curve.
    ///
    /// This constructor is intentionally low-level: it does not try to infer
    /// special forms such as constants or pure `y`-multiples, it simply stores
    /// the caller-supplied pair in the standard decomposition
    /// `F(E) = F(x) ⊕ yF(x)`.
    pub fn new(
        curve: ShortWeierstrassCurve<F>,
        a_part: RationalFunction<F>,
        b_part: RationalFunction<F>,
    ) -> Self {
        Self {
            curve,
            a_part,
            b_part,
        }
    }

    /// Embeds a rational function `A(x)` as `A(x) + y * 0`.
    pub fn from_rational_function(
        curve: ShortWeierstrassCurve<F>,
        function: RationalFunction<F>,
    ) -> Self {
        Self::new(curve, function, RationalFunction::constant(F::zero()))
    }

    pub fn zero(curve: ShortWeierstrassCurve<F>) -> Self {
        Self::constant(curve, F::zero())
    }

    pub fn one(curve: ShortWeierstrassCurve<F>) -> Self {
        Self::constant(curve, F::one())
    }

    pub fn y(curve: ShortWeierstrassCurve<F>) -> Self {
        Self::new(
            curve,
            RationalFunction::constant(F::zero()),
            RationalFunction::constant(F::one()),
        )
    }

    pub(crate) fn constant(curve: ShortWeierstrassCurve<F>, value: F::Elem) -> Self {
        Self::from_rational_function(curve, RationalFunction::constant(value))
    }

    pub(crate) fn from_finite_affine_point(
        curve: ShortWeierstrassCurve<F>,
        point: &AffinePoint<F>,
        use_x: bool,
    ) -> Self {
        let coordinate = match point {
            AffinePoint::Infinity => {
                panic!("finite-point coordinate embedding requires a finite affine point")
            }
            AffinePoint::Finite { x, y } => {
                if use_x {
                    x.clone()
                } else {
                    y.clone()
                }
            }
        };

        Self::constant(curve, coordinate)
    }
}
