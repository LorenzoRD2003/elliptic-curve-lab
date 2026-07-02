use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve, short_weierstrass::function_fields::ShortWeierstrassFunction,
};
use crate::fields::traits::*;

/// Point of a short-Weierstrass curve with coordinates in its own function
/// field `F(E)`.
///
/// This is the function-field analogue of [`crate::elliptic_curves::AffinePoint`]:
/// it stores either the point at infinity or one affine pair of function-field coordinates.
///
/// The main educational use of this type is to represent the generic point
/// `P_gen = (x, y)` and its images under rational group-law formulas such as
/// translation by a fixed point, doubling, or later scalar multiplication.
#[derive(Debug, Clone)]
pub enum ShortWeierstrassFunctionFieldPoint<F: Field> {
    /// The point at infinity.
    Infinity,
    /// An affine point with coordinates in `F(E)`.
    Affine {
        x: ShortWeierstrassFunction<F>,
        y: ShortWeierstrassFunction<F>,
    },
}

impl<F: Field> ShortWeierstrassFunctionFieldPoint<F> {
    /// Returns the point at infinity.
    pub fn infinity() -> Self {
        Self::Infinity
    }

    /// Builds an affine function-field point after checking that both
    /// coordinates live on the same ambient curve and satisfy its equation.
    pub fn affine(
        x: ShortWeierstrassFunction<F>,
        y: ShortWeierstrassFunction<F>,
    ) -> Result<Self, CurveError> {
        Self::validate_affine_coordinates_on_curve(&x, &y, x.curve())?;

        Ok(Self::Affine { x, y })
    }

    /// Returns whether this point is the point at infinity.
    pub fn is_infinity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Returns the affine `x`-coordinate when it exists.
    pub fn x(&self) -> Option<&ShortWeierstrassFunction<F>> {
        match self {
            Self::Infinity => None,
            Self::Affine { x, .. } => Some(x),
        }
    }

    /// Returns the affine `y`-coordinate when it exists.
    pub fn y(&self) -> Option<&ShortWeierstrassFunction<F>> {
        match self {
            Self::Infinity => None,
            Self::Affine { y, .. } => Some(y),
        }
    }

    /// Returns the underlying curve when the point is affine.
    pub fn curve(&self) -> Option<&ShortWeierstrassCurve<F>> {
        self.x().map(ShortWeierstrassFunction::curve)
    }
}

impl<F: Field> PartialEq for ShortWeierstrassFunctionFieldPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Affine {
                    x: left_x,
                    y: left_y,
                },
                Self::Affine {
                    x: right_x,
                    y: right_y,
                },
            ) => left_x == right_x && left_y == right_y,
            _ => false,
        }
    }
}
