use crate::elliptic_curves::CurveError;
use crate::fields::SqrtField;

use super::CurveModel;

/// Curve models that admit affine coordinate validation.
pub trait AffineCurveModel: CurveModel {
    /// Builds a point from affine coordinates after checking that it lies on
    /// the curve.
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError>;
}

/// Curve models that can lift an `x`-coordinate into one or two affine points.
///
/// This trait is intentionally about the *shape* of the curve equation rather
/// than about enumeration or group-law operations. It models the common
/// situation where the curve equation determines `y^2` from a chosen `x`.
pub trait LiftXCoordinate: AffineCurveModel
where
    Self::BaseField: SqrtField<Elem = Self::Elem>,
{
    /// Returns the right-hand side of the curve equation as a function of `x`.
    fn rhs(&self, x: &Self::Elem) -> Self::Elem;

    /// Builds one point above the given `x` when a square root exists.
    ///
    /// Which square root is chosen is delegated to the base field's
    /// [`SqrtField`] implementation.
    fn point_from_x(&self, x: Self::Elem) -> Option<Self::Point> {
        let y = Self::BaseField::sqrt(&self.rhs(&x))?;
        self.point(x, y).ok()
    }

    /// Builds the two points above the given `x` when square roots exist.
    ///
    /// When the only root is `0`, both returned points are the same because
    /// the two square roots coincide.
    fn points_from_x(&self, x: Self::Elem) -> Option<(Self::Point, Self::Point)> {
        let (left_y, right_y) = Self::BaseField::sqrt_pair(&self.rhs(&x))?;
        let left = self.point(x.clone(), left_y).ok()?;
        let right = self.point(x, right_y).ok()?;
        Some((left, right))
    }
}
