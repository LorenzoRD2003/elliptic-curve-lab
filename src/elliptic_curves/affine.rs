use core::fmt;
use std::hash::{Hash, Hasher};

use crate::fields::Field;

/// Affine point representation for elliptic curves.
///
/// The point at infinity is modeled explicitly instead of storing placeholder
/// affine coordinates alongside a boolean flag.
pub enum AffinePoint<F: Field> {
    Infinity,
    Finite { x: F::Elem, y: F::Elem },
}

impl<F: Field> AffinePoint<F> {
    /// Builds a finite affine point.
    pub fn new(x: F::Elem, y: F::Elem) -> Self {
        Self::Finite { x, y }
    }

    /// Builds the distinguished point at infinity.
    pub fn infinity() -> Self {
        Self::Infinity
    }

    /// Returns a compact affine-coordinate string or `O` for the identity.
    pub fn to_coordinates_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        match self {
            Self::Infinity => "O".to_string(),
            Self::Finite { x, y } => format!("({x}, {y})"),
        }
    }

    /// Returns the affine negation `(x, -y)` and fixes the identity.
    ///
    /// This is the standard coordinate involution used by curve models whose
    /// defining equation is symmetric under `y -> -y`, such as short
    /// Weierstrass form.
    pub fn neg(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Finite { x, y } => Self::Finite {
                x: x.clone(),
                y: F::neg(y),
            },
        }
    }

    /// Transports the stored affine coordinates through a caller-supplied map.
    ///
    /// This helper acts only on the point representation:
    ///
    /// - `O` stays `O`
    /// - a finite point `(x, y)` is sent to `(f(x), f(y))`
    ///
    /// It does not certify that the image lies on any particular target curve.
    /// Callers that care about curve membership or compatibility with a target
    /// model must validate that separately.
    pub fn map_coordinates<G, M>(&self, mut map: M) -> AffinePoint<G>
    where
        G: Field,
        M: FnMut(&F::Elem) -> G::Elem,
    {
        match self {
            Self::Infinity => AffinePoint::Infinity,
            Self::Finite { x, y } => AffinePoint::new(map(x), map(y)),
        }
    }

    /// Returns whether this point is the distinguished identity element.
    pub fn is_identity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Returns affine coordinate references when the point is finite.
    pub(crate) fn finite_coordinates(point: &Self) -> Option<(&F::Elem, &F::Elem)> {
        match point {
            Self::Infinity => None,
            Self::Finite { x, y } => Some((x, y)),
        }
    }

    /// Returns the `x`-coordinate when the point is finite.
    pub(crate) fn x_coordinate(point: &Self) -> Option<&F::Elem> {
        Self::finite_coordinates(point).map(|(x, _)| x)
    }

    /// Returns the `y`-coordinate when the point is finite.
    pub(crate) fn y_coordinate(point: &Self) -> Option<&F::Elem> {
        Self::finite_coordinates(point).map(|(_, y)| y)
    }
}

impl<F: Field> PartialEq for AffinePoint<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Finite {
                    x: x_left,
                    y: y_left,
                },
                Self::Finite {
                    x: x_right,
                    y: y_right,
                },
            ) => F::eq(x_left, x_right) && F::eq(y_left, y_right),
            _ => false,
        }
    }
}

impl<F: Field> Eq for AffinePoint<F> {}

impl<F> Hash for AffinePoint<F>
where
    F: Field,
    F::Elem: Hash,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        match self {
            Self::Infinity => {
                0_u8.hash(state);
            }
            Self::Finite { x, y } => {
                1_u8.hash(state);
                x.hash(state);
                y.hash(state);
            }
        }
    }
}

impl<F: Field> Clone for AffinePoint<F> {
    fn clone(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Finite { x, y } => Self::Finite {
                x: x.clone(),
                y: y.clone(),
            },
        }
    }
}

impl<F> fmt::Display for AffinePoint<F>
where
    F: Field,
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_coordinates_string())
    }
}

impl<F: Field> fmt::Debug for AffinePoint<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Infinity => write!(f, "AffinePoint::Infinity"),
            Self::Finite { x, y } => f
                .debug_struct("AffinePoint")
                .field("x", x)
                .field("y", y)
                .finish(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::AffinePoint;
    use crate::fields::{Field, Fp};

    type F7 = Fp<7>;

    #[test]
    fn finite_constructor_marks_point_as_finite() {
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));

        match point {
            AffinePoint::Finite { x, y } => {
                assert!(F7::eq(&x, &F7::from_i64(2)));
                assert!(F7::eq(&y, &F7::from_i64(5)));
            }
            AffinePoint::Infinity => panic!("expected a finite affine point"),
        }
    }

    #[test]
    fn infinity_constructor_builds_distinguished_variant() {
        let point = AffinePoint::<F7>::infinity();

        assert!(matches!(point, AffinePoint::Infinity));
        assert!(point.is_identity());
    }

    #[test]
    fn finite_points_are_not_identity() {
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));

        assert!(!point.is_identity());
    }

    #[test]
    fn negates_finite_points_by_flipping_the_y_coordinate() {
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));
        let negated = point.neg();

        match negated {
            AffinePoint::Finite { x, y } => {
                assert!(F7::eq(&x, &F7::from_i64(2)));
                assert!(F7::eq(&y, &F7::from_i64(-5)));
            }
            AffinePoint::Infinity => panic!("expected a finite affine point"),
        }
    }

    #[test]
    fn negation_fixes_the_identity() {
        let point = AffinePoint::<F7>::infinity();

        assert_eq!(point.neg(), AffinePoint::Infinity);
    }

    #[test]
    fn coordinate_map_preserves_the_identity() {
        let point = AffinePoint::<F7>::infinity();
        let image = point.map_coordinates::<F7, _>(|coordinate| *coordinate);

        assert_eq!(image, AffinePoint::Infinity);
    }

    #[test]
    fn coordinate_map_transforms_finite_coordinates_componentwise() {
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));
        let image = point.map_coordinates::<F7, _>(|coordinate| F7::add(coordinate, &F7::one()));

        assert_eq!(
            image,
            AffinePoint::<F7>::new(F7::from_i64(3), F7::from_i64(6))
        );
    }

    #[test]
    fn coordinate_helpers_return_none_for_identity() {
        let point = AffinePoint::<F7>::infinity();

        assert_eq!(AffinePoint::<F7>::finite_coordinates(&point), None);
        assert_eq!(AffinePoint::<F7>::x_coordinate(&point), None);
        assert_eq!(AffinePoint::<F7>::y_coordinate(&point), None);
    }

    #[test]
    fn coordinate_helpers_expose_finite_coordinates() {
        let point = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(5));

        let (x, y) = AffinePoint::<F7>::finite_coordinates(&point).expect("point should be finite");
        assert!(F7::eq(x, &F7::from_i64(2)));
        assert!(F7::eq(y, &F7::from_i64(5)));
        assert!(F7::eq(
            AffinePoint::<F7>::x_coordinate(&point).expect("x should exist"),
            &F7::from_i64(2)
        ));
        assert!(F7::eq(
            AffinePoint::<F7>::y_coordinate(&point).expect("y should exist"),
            &F7::from_i64(5)
        ));
    }
}
