use crate::elliptic_curves::{AffinePoint, TwistedEdwardsCurve};
use crate::fields::traits::*;

impl<F: Field> TwistedEdwardsCurve<F> {
    pub(crate) fn contains_affine_point(&self, point: &AffinePoint<F>) -> bool {
        match point {
            AffinePoint::Infinity => false,
            AffinePoint::Finite { x, y } => {
                let x_sq = F::square(x);
                let y_sq = F::square(y);
                let left = F::add(&F::mul(self.a(), &x_sq), &y_sq);
                let right = F::add(&F::one(), &F::mul(self.d(), &F::mul(&x_sq, &y_sq)));
                F::eq(&left, &right)
            }
        }
    }
}
