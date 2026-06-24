use core::fmt;

use crate::elliptic_curves::{AffinePoint, CurveError};
use crate::fields::traits::Field;

/// Projective `x`-line coordinate for Montgomery ladder work.
///
/// Finite values are stored as projective pairs `(X : Z)` and interpreted as
/// the affine ratio `x = X / Z` when `Z ≠ 0`.
///
/// The distinguished point at infinity is modeled explicitly instead of being
/// encoded by a degenerate `Z = 0` pair.
pub enum MontgomeryXzPoint<F: Field> {
    Infinity,
    Finite { x: F::Elem, z: F::Elem },
}

impl<F: Field> MontgomeryXzPoint<F> {
    /// Builds one finite `X:Z` representative without checking whether `Z` is
    /// invertible.
    pub fn new(x: F::Elem, z: F::Elem) -> Self {
        Self::Finite { x, z }
    }

    /// Builds the distinguished point at infinity on the `x`-line.
    pub fn infinity() -> Self {
        Self::Infinity
    }

    /// Lifts one affine `x` value into the normalized projective chart `Z = 1`.
    pub fn from_affine_x(x: F::Elem) -> Self {
        Self::Finite { x, z: F::one() }
    }

    /// Forgets the sign of `y` and keeps only the `x`-line class of one affine
    /// Montgomery point.
    pub fn from_affine_point(point: &AffinePoint<F>) -> Self {
        match point {
            AffinePoint::Infinity => Self::Infinity,
            AffinePoint::Finite { x, .. } => Self::from_affine_x(x.clone()),
        }
    }

    /// Builds one finite `X:Z` representative unless `Z = 0`, in which case
    /// the result is treated as the explicit point at infinity on the `x`-line.
    pub fn from_xz_or_infinity(x: F::Elem, z: F::Elem) -> Self {
        if F::is_zero(&z) {
            Self::Infinity
        } else {
            Self::new(x, z)
        }
    }

    /// Returns whether this value is the explicit point at infinity.
    pub fn is_infinity(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    /// Returns whether the stored finite representative already lies in the
    /// normalized `Z = 1` chart.
    ///
    /// The point at infinity is treated as already normalized.
    pub fn is_normalized(&self) -> bool {
        match self {
            Self::Infinity => true,
            Self::Finite { z, .. } => F::eq(z, &F::one()),
        }
    }

    /// Returns borrowed finite coordinates when this value stores one finite
    /// `X:Z` representative.
    pub fn finite_coordinates(&self) -> Option<(&F::Elem, &F::Elem)> {
        match self {
            Self::Infinity => None,
            Self::Finite { x, z } => Some((x, z)),
        }
    }

    /// Recovers the affine `x` value.
    ///
    /// Finite representatives require `Z` to be invertible. The point at
    /// infinity returns `None`.
    pub fn to_affine_x(&self) -> Result<Option<F::Elem>, CurveError> {
        match self {
            Self::Infinity => Ok(None),
            Self::Finite { x, z } => {
                let z_inverse = F::inverse(z).map_err(CurveError::Field)?;
                Ok(Some(F::mul(x, &z_inverse)))
            }
        }
    }

    /// Returns an equivalent representative with `Z = 1`.
    pub fn normalize(&self) -> Result<Self, CurveError> {
        match self {
            Self::Infinity => Ok(Self::Infinity),
            Self::Finite { x, z } => {
                let z_inverse = F::inverse(z).map_err(CurveError::Field)?;
                Ok(Self::Finite {
                    x: F::mul(x, &z_inverse),
                    z: F::one(),
                })
            }
        }
    }

    /// Returns whether both values store exactly the same representative.
    pub fn has_same_representative_as(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Finite {
                    x: left_x,
                    z: left_z,
                },
                Self::Finite {
                    x: right_x,
                    z: right_z,
                },
            ) => F::eq(left_x, right_x) && F::eq(left_z, right_z),
            _ => false,
        }
    }

    /// Returns a compact coordinate string or `O_x` for infinity.
    pub fn to_coordinates_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        match self {
            Self::Infinity => "O_x".to_string(),
            Self::Finite { x, z } => format!("({x} : {z})"),
        }
    }
}

impl<F: Field> Clone for MontgomeryXzPoint<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        match self {
            Self::Infinity => Self::Infinity,
            Self::Finite { x, z } => Self::Finite {
                x: x.clone(),
                z: z.clone(),
            },
        }
    }
}

impl<F: Field> PartialEq for MontgomeryXzPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Infinity, Self::Infinity) => true,
            (
                Self::Finite {
                    x: left_x,
                    z: left_z,
                },
                Self::Finite {
                    x: right_x,
                    z: right_z,
                },
            ) => {
                if F::is_zero(left_z) || F::is_zero(right_z) {
                    return self.has_same_representative_as(other);
                }

                F::eq(&F::mul(left_x, right_z), &F::mul(right_x, left_z))
            }
            _ => false,
        }
    }
}

impl<F: Field> Eq for MontgomeryXzPoint<F> {}

impl<F: Field> fmt::Display for MontgomeryXzPoint<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.to_coordinates_string())
    }
}

impl<F: Field> fmt::Debug for MontgomeryXzPoint<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Infinity => formatter.write_str("MontgomeryXzPoint::Infinity"),
            Self::Finite { x, z } => formatter
                .debug_struct("MontgomeryXzPoint")
                .field("x", x)
                .field("z", z)
                .finish(),
        }
    }
}
