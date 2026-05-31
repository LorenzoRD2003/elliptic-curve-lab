use crate::fields::Field;

/// Affine point representation for elliptic curves.
///
/// The point at infinity is modeled explicitly instead of storing placeholder
/// affine coordinates alongside a boolean flag.
#[derive(Clone, Debug)]
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

    /// Returns whether this point is the distinguished identity element.
    pub fn is_identity(&self) -> bool {
        matches!(self, Self::Infinity)
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

#[cfg(test)]
mod tests {
    use super::AffinePoint;
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
}
