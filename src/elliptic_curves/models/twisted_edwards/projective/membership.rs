use crate::elliptic_curves::{
    TwistedEdwardsCurve, twisted_edwards::projective::ExtendedTwistedEdwardsPoint,
};
use crate::fields::traits::*;

impl<F: Field> TwistedEdwardsCurve<F> {
    /// Returns whether one extended point satisfies the twisted-Edwards
    /// projective equations
    ///
    /// `aX^2 + Y^2 = Z^2 + dT^2` and `XY = ZT`.
    ///
    /// The all-zero tuple is rejected because it does not define a projective
    /// point. Points with `Z = 0` are allowed here: they belong to the
    /// projective closure even though they do not admit affine recovery.
    pub(crate) fn contains_extended_point(&self, point: &ExtendedTwistedEdwardsPoint<F>) -> bool {
        if point.is_zero_tuple() {
            return false;
        }
        let left = F::add(
            &F::mul(self.a(), &F::square(point.x())),
            &F::square(point.y()),
        );
        let right = F::add(
            &F::square(point.z()),
            &F::mul(self.d(), &F::square(point.t())),
        );
        let structural = F::eq(&F::mul(point.x(), point.y()), &F::mul(point.z(), point.t()));
        structural && F::eq(&left, &right)
    }
}
