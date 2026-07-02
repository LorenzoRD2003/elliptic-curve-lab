use crate::elliptic_curves::{AffinePoint, MontgomeryCurve};
use crate::fields::traits::*;

impl<F: Field> MontgomeryCurve<F> {
    pub(crate) fn contains_affine_point(&self, point: &AffinePoint<F>) -> bool {
        match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => {
                let left = F::mul(self.b(), &F::square(y));
                let right = self.rhs_value(x);
                F::eq(&left, &right)
            }
        }
    }
}
