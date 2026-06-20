use crate::elliptic_curves::{CurveError, traits::CurveModel};

/// Curve models that admit affine coordinate validation.
pub trait AffineCurveModel: CurveModel {
    /// Builds a point from affine coordinates after checking that it lies on
    /// the curve.
    fn point(&self, x: Self::Elem, y: Self::Elem) -> Result<Self::Point, CurveError>;
}

/// The finite fiber of the projection `x : E -> A^1` over one chosen base-field
/// `x`-coordinate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LiftedPoints<P> {
    NoPoint,
    OnePoint(P),
    TwoPoints(P, P),
}

/// Curve models that can lift an `x`-coordinate into its affine fiber.
///
/// This trait is intentionally about the *shape* of the curve equation rather
/// than about enumeration or group-law operations. It models the fiber of the
/// coordinate projection over one chosen `x`, whether that fiber is recovered
/// from a square root, a shifted quadratic solve, or another model-specific
/// route.
pub trait LiftXCoordinate: AffineCurveModel {
    /// Returns the affine fiber above the chosen `x`.
    fn lift_x(&self, x: Self::Elem) -> Result<LiftedPoints<Self::Point>, CurveError>;

    /// Returns one point above the chosen `x` when the fiber is non-empty.
    fn point_from_x(&self, x: Self::Elem) -> Result<Option<Self::Point>, CurveError> {
        match self.lift_x(x)? {
            LiftedPoints::NoPoint => Ok(None),
            LiftedPoints::OnePoint(point) | LiftedPoints::TwoPoints(point, _) => Ok(Some(point)),
        }
    }

    /// Returns the fiber above the chosen `x` as an optional pair.
    ///
    /// This keeps the old ergonomic surface for callers that naturally expect
    /// either “no point” or “one/two points”, while the canonical trait method
    /// remains [`Self::lift_x`].
    fn points_from_x(
        &self,
        x: Self::Elem,
    ) -> Result<Option<(Self::Point, Self::Point)>, CurveError> {
        match self.lift_x(x)? {
            LiftedPoints::NoPoint => Ok(None),
            LiftedPoints::OnePoint(point) => Ok(Some((point.clone(), point))),
            LiftedPoints::TwoPoints(left, right) => Ok(Some((left, right))),
        }
    }
}
