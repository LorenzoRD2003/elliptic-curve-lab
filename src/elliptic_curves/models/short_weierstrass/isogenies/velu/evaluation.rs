use std::hash::Hash;

use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::isogenies::{VeluIsogeny, velu::VeluKernelData},
    traits::{AffineCurveModel, GroupCurveModel},
};
use crate::fields::traits::Field;
use crate::isogenies::error::IsogenyError;

type TranslationSumCoordinates<F> = Option<(<F as Field>::Elem, <F as Field>::Elem)>;

impl<F: Field> VeluKernelData<F>
where
    F::Elem: Clone + Eq + Hash,
{
    pub(super) fn translation_correction_sums(
        &self,
        domain: &ShortWeierstrassCurve<F>,
        point: &AffinePoint<F>,
    ) -> Result<(F::Elem, F::Elem), IsogenyError> {
        debug_assert!(
            AffinePoint::finite_coordinates(point).is_some(),
            "finite non-kernel translation sums should expose a finite point"
        );
        let x = AffinePoint::x_coordinate(point)
            .expect("finite non-kernel translation sums should expose an affine x-coordinate");
        let y = AffinePoint::y_coordinate(point)
            .expect("finite non-kernel translation sums should expose an affine y-coordinate");
        let mut x_sum = x.clone();
        let mut y_sum = y.clone();

        for kernel_point in self.kernel_nonzero_points() {
            let translated = domain.add(point, kernel_point)?;

            let translated_coordinates = AffinePoint::finite_coordinates(&translated);
            let (translated_x, translated_y) = translated_coordinates.expect(
                "if P is outside the kernel, then P + Q cannot be the identity for non-zero Q in the kernel",
            );

            let kernel_coordinates = AffinePoint::finite_coordinates(kernel_point);
            let (kernel_x, kernel_y) =
                kernel_coordinates.expect("kernel_nonzero_points never contains the identity");

            x_sum = F::add(&x_sum, &F::sub(translated_x, kernel_x));
            y_sum = F::add(&y_sum, &F::sub(translated_y, kernel_y));
        }

        Ok((x_sum, y_sum))
    }
}

impl<F: Field> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F::Elem: Clone + Eq + Hash,
{
    pub(crate) fn evaluate_non_kernel_point(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, IsogenyError> {
        let finite_point = self
            .require_non_kernel_finite_point(point)?
            .expect("non-kernel short-Weierstrass Vélu evaluation should be finite");
        let (x, y) = VeluKernelData::from_kernel(&self.domain, &self.kernel)
            .translation_correction_sums(&self.domain, finite_point)?;
        self.codomain.point(x, y).map_err(IsogenyError::from)
    }

    /// Computes the affine-coordinate translation sums for Vélu's map.
    ///
    /// If `G* = G \ {O}` denotes the non-identity kernel points, the intended
    /// implementation of Vélu uses the classical affine formulas
    ///
    /// `x_φ(P) = x(P) + Σ_{Q in G*} (x(P + Q) - x(Q))`
    /// `y_φ(P) = y(P) + Σ_{Q in G*} (y(P + Q) - y(Q))`.
    ///
    /// This helper computes exactly those sums for a finite affine point `P`
    /// outside the kernel. If `P` lies in the kernel, then the Vélu image is
    /// the identity in the codomain, so there are no affine output
    /// coordinates to return and the method yields `Ok(None)`.
    pub(crate) fn translation_sum_coordinates(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<TranslationSumCoordinates<F>, IsogenyError> {
        match self.require_non_kernel_finite_point(point)? {
            Some(finite_point) => VeluKernelData::from_kernel(&self.domain, &self.kernel)
                .translation_correction_sums(&self.domain, finite_point)
                .map(Some),
            None => Ok(None),
        }
    }
}
