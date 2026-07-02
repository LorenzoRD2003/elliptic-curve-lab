use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::{AffinePoint, CurveError};

/// Extended twisted-Edwards point `(X : Y : Z : T)`.
///
/// This representation is attached to the affine twisted-Edwards chart
/// `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2` through
///
/// `x = X / Z`, `y = Y / Z`, `T = XY / Z`.
///
/// The extra coordinate is constrained by the structural identity `XY = ZT`.
///
/// For curve points the homogeneous twisted-Edwards equation becomes
/// `aX^2 + Y^2 = Z^2 + dT^2`.
///
/// Affine points embed into the normalized form `Z = 1` by
/// `(x, y) -> (x : y : 1 : xy)`, where the neutral element is
/// `(0 : 1 : 1 : 0)`.
pub struct ExtendedTwistedEdwardsPoint<F: Field> {
    x: F::Elem,
    y: F::Elem,
    z: F::Elem,
    t: F::Elem,
}

impl<F: Field> ExtendedTwistedEdwardsPoint<F> {
    /// Builds one raw extended-coordinate representative without validating
    /// curve membership.
    pub fn new(x: F::Elem, y: F::Elem, z: F::Elem, t: F::Elem) -> Self {
        Self { x, y, z, t }
    }

    /// Returns the neutral element `(0 : 1 : 1 : 0)`.
    pub fn identity() -> Self {
        Self::new(F::zero(), F::one(), F::one(), F::zero())
    }

    /// Embeds one finite affine point into the normalized `Z = 1` chart.
    ///
    /// Twisted Edwards has finite affine identity, so `AffinePoint::Infinity`
    /// is rejected honestly here.
    pub fn from_affine(point: &AffinePoint<F>) -> Result<Self, CurveError> {
        match point {
            AffinePoint::Infinity => Err(CurveError::PointNotOnCurve),
            AffinePoint::Finite { x, y } => {
                Ok(Self::new(x.clone(), y.clone(), F::one(), F::mul(x, y)))
            }
        }
    }

    /// Recovers the affine coordinates `x = X/Z`, `y = Y/Z`.
    ///
    /// This requires `Z` to be invertible. If `Z = 0`, the method returns the
    /// underlying field inversion error instead of inventing an affine fallback.
    pub fn to_affine(&self) -> Result<AffinePoint<F>, CurveError> {
        let z_inverse = F::inverse(&self.z).map_err(CurveError::Field)?;
        Ok(AffinePoint::new(
            F::mul(&self.x, &z_inverse),
            F::mul(&self.y, &z_inverse),
        ))
    }

    /// Returns the additive inverse `(-X : Y : Z : -T)`.
    ///
    /// This is the projective lift of the affine twisted-Edwards involution
    /// `-(x, y) = (-x, y)`.
    pub fn neg(&self) -> Self {
        Self::new(
            F::neg(&self.x),
            self.y.clone(),
            self.z.clone(),
            F::neg(&self.t),
        )
    }

    /// Returns whether this representative is projectively equal to the
    /// neutral element `(0 : 1 : 1 : 0)`.
    pub fn is_identity(&self) -> bool {
        self == &Self::identity()
    }

    /// Returns the stored `X` coordinate.
    pub fn x(&self) -> &F::Elem {
        &self.x
    }

    /// Returns the stored `Y` coordinate.
    pub fn y(&self) -> &F::Elem {
        &self.y
    }

    /// Returns the stored `Z` coordinate.
    pub fn z(&self) -> &F::Elem {
        &self.z
    }

    /// Returns the stored `T` coordinate.
    pub fn t(&self) -> &F::Elem {
        &self.t
    }

    /// Returns whether both points store exactly the same representative.
    pub fn has_same_representative_as(&self, other: &Self) -> bool {
        F::eq(&self.x, &other.x)
            && F::eq(&self.y, &other.y)
            && F::eq(&self.z, &other.z)
            && F::eq(&self.t, &other.t)
    }

    pub(crate) fn is_zero_tuple(&self) -> bool {
        F::is_zero(&self.x) && F::is_zero(&self.y) && F::is_zero(&self.z) && F::is_zero(&self.t)
    }
}

impl<F: Field> Clone for ExtendedTwistedEdwardsPoint<F> {
    fn clone(&self) -> Self {
        Self::new(
            self.x.clone(),
            self.y.clone(),
            self.z.clone(),
            self.t.clone(),
        )
    }
}

impl<F: Field> PartialEq for ExtendedTwistedEdwardsPoint<F> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero_tuple() || other.is_zero_tuple() {
            return self.has_same_representative_as(other);
        }

        F::eq(&F::mul(&self.x, &other.y), &F::mul(&other.x, &self.y))
            && F::eq(&F::mul(&self.x, &other.z), &F::mul(&other.x, &self.z))
            && F::eq(&F::mul(&self.x, &other.t), &F::mul(&other.x, &self.t))
            && F::eq(&F::mul(&self.y, &other.z), &F::mul(&other.y, &self.z))
            && F::eq(&F::mul(&self.y, &other.t), &F::mul(&other.y, &self.t))
            && F::eq(&F::mul(&self.z, &other.t), &F::mul(&other.z, &self.t))
    }
}

impl<F: Field> Eq for ExtendedTwistedEdwardsPoint<F> {}

impl<F: Field> fmt::Display for ExtendedTwistedEdwardsPoint<F>
where
    F::Elem: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} : {} : {} : {})", self.x, self.y, self.z, self.t)
    }
}

impl<F: Field> fmt::Debug for ExtendedTwistedEdwardsPoint<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ExtendedTwistedEdwardsPoint")
            .field("x", &self.x)
            .field("y", &self.y)
            .field("z", &self.z)
            .field("t", &self.t)
            .finish()
    }
}
