use std::hash::Hash;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::fields::Field;

use crate::isogenies::{IsogenyError, IsogenyKernel};

use crate::isogenies::velu::short_weierstrass::VeluKernelData;

impl<F> VeluKernelData<F>
where
    F: Field,
    F::Elem: Clone + Eq + Hash,
{
    pub(super) fn from_kernel(
        domain: &ShortWeierstrassCurve<F>,
        kernel: &IsogenyKernel<ShortWeierstrassCurve<F>>,
    ) -> Self {
        let a = domain.a().clone();
        let b = domain.b().clone();

        let kernel_nonzero_points = kernel
            .points()
            .split_first()
            .map(|(_, tail)| tail.to_vec())
            .unwrap_or_default();
        let mut a_correction_sum = F::zero();
        let mut b_correction_sum = F::zero();

        for kernel_point in &kernel_nonzero_points {
            let kernel_x = AffinePoint::x_coordinate(kernel_point);
            debug_assert!(
                kernel_x.is_some(),
                "kernel_nonzero_points never contains the identity"
            );
            let x = kernel_x.expect("kernel_nonzero_points never contains the identity");

            let x_squared = F::square(x);
            let x_cubed = F::cube(x);

            let a_term = F::add(&F::mul(&F::from_i64(3), &x_squared), &a);
            a_correction_sum = F::add(&a_correction_sum, &a_term);

            let five_x_cubed = F::mul(&F::from_i64(5), &x_cubed);
            let three_a_x = F::mul(&F::from_i64(3), &F::mul(&a, x));
            let two_b = F::mul(&F::from_i64(2), &b);
            let b_term = F::add(&F::add(&five_x_cubed, &three_a_x), &two_b);
            b_correction_sum = F::add(&b_correction_sum, &b_term);
        }

        Self {
            kernel_nonzero_points,
            a_correction_sum,
            b_correction_sum,
        }
    }

    /// Computes the short-Weierstrass codomain curve from Vélu's coefficient formulas.
    ///
    /// For a finite kernel `G` with `G* = G \ {O}`, the coefficient updates are
    ///
    /// `a' = a - 5 * sum_{Q in G*} (3x_Q^2 + a)`
    /// `b' = b - 7 * sum_{Q in G*} (5x_Q^3 + 3ax_Q + 2b)`.
    pub(super) fn codomain_curve(
        &self,
        domain: &ShortWeierstrassCurve<F>,
    ) -> Result<ShortWeierstrassCurve<F>, IsogenyError> {
        let a = domain.a().clone();
        let b = domain.b().clone();
        let a_prime = F::sub(&a, &F::mul(&F::from_i64(5), &self.a_correction_sum));
        let b_prime = F::sub(&b, &F::mul(&F::from_i64(7), &self.b_correction_sum));

        ShortWeierstrassCurve::new(a_prime, b_prime).map_err(IsogenyError::from)
    }
}
