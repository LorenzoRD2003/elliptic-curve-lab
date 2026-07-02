use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve};
use crate::fields::traits::*;

impl<F: Field> GeneralWeierstrassCurve<F> {
    pub(crate) fn contains_affine_point(&self, point: &AffinePoint<F>) -> bool {
        match point {
            AffinePoint::Infinity => true,
            AffinePoint::Finite { x, y } => {
                let left = F::add(
                    &F::add(&F::square(y), &F::mul(&F::mul(self.a1(), x), y)),
                    &F::mul(self.a3(), y),
                );
                let right = F::add(
                    &F::add(&F::cube(x), &F::mul(self.a2(), &F::square(x))),
                    &F::add(&F::mul(self.a4(), x), self.a6()),
                );
                F::eq(&left, &right)
            }
        }
    }
}
