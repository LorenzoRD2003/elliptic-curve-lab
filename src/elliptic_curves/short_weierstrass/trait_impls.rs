use core::fmt;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::frobenius::relative_frobenius_point;
use crate::elliptic_curves::invariants::HasJInvariant;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModel, LiftXCoordinate, RelativeFrobeniusCurveModel,
};
use crate::fields::{Field, SqrtField};

impl<F: Field> fmt::Display for ShortWeierstrassCurve<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_equation_string())
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassCurve<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShortWeierstrassCurve")
            .field(
                "equation",
                &format_args!("y^2 = x^3 + ({:?})x + ({:?})", self.a, self.b),
            )
            .field("a", &self.a)
            .field("b", &self.b)
            .finish()
    }
}

impl<F: Field> CurveModel for ShortWeierstrassCurve<F> {
    type Elem = F::Elem;
    type BaseField = F;
    type Point = AffinePoint<F>;

    fn identity(&self) -> Self::Point {
        AffinePoint::infinity()
    }

    fn is_identity(&self, point: &Self::Point) -> bool {
        point.is_identity()
    }

    fn contains(&self, point: &Self::Point) -> bool {
        match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => {
                let left = F::square(y);
                let right = self.rhs_value(x);
                F::eq(&left, &right)
            }
        }
    }
}

impl<F: Field> HasJInvariant for ShortWeierstrassCurve<F> {
    fn j_invariant(&self) -> Self::Elem {
        ShortWeierstrassCurve::j_invariant(self)
    }
}

impl<F: Field> AffineCurveModel for ShortWeierstrassCurve<F> {
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError> {
        let point = self.unchecked_point(x, y);
        if self.contains(&point) {
            Ok(point)
        } else {
            Err(CurveError::PointNotOnCurve)
        }
    }
}

impl<F: SqrtField> LiftXCoordinate for ShortWeierstrassCurve<F> {
    fn rhs(&self, x: &Self::Elem) -> Self::Elem {
        self.rhs_value(x)
    }
}

impl<F: crate::fields::FiniteField> RelativeFrobeniusCurveModel for ShortWeierstrassCurve<F> {
    fn relative_frobenius(&self, point: &Self::Point) -> Result<Self::Point, CurveError> {
        relative_frobenius_point(self, point)
    }
}
