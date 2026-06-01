mod codomain;
mod dual;
mod evaluation;

use std::hash::Hash;

use crate::elliptic_curves::{AffineCurveModel, AffinePoint, CurveModel, ShortWeierstrassCurve};
use crate::fields::Field;

use crate::isogenies::{IsogenyError, IsogenyKernel, VeluIsogeny};

use super::SupportsVeluConstruction;

pub use dual::{DualVeluIsogeny, verify_left_dual_relation, verify_right_dual_relation};

/// Internal precomputation bucket for the short-Weierstrass Vélu construction.
///
/// Mathematically, Vélu's formulas are driven by one finite subgroup
/// `G \subset E` and, more precisely, by the non-identity subset
/// `G* = G \ {O}`. Both parts of the construction that matter here
///
/// - the codomain coefficients `a'` and `b'`
/// - the affine translation sums used to evaluate points outside the kernel
///
/// are built from repeated passes over those same kernel points.
///
/// This struct exists to keep that dependency explicit. Instead of letting the
/// codomain logic and the evaluation logic each derive their own ad hoc
/// intermediate values from the kernel, `VeluKernelData` stores the shared
/// finite data extracted from `G*` once and then lets both sides consume the
/// same source of truth.
///
/// In particular:
///
/// - `kernel_nonzero_points` records the list of non-identity kernel points
///   that appear in the classical Vélu sums;
/// - `a_correction_sum` stores the accumulated
///   `sum_{Q in G*} (3x_Q^2 + a)` term used in the `a'` update;
/// - `b_correction_sum` stores the accumulated
///   `sum_{Q in G*} (5x_Q^3 + 3ax_Q + 2b)` term used in the `b'` update.
#[derive(Clone, Debug)]
struct VeluKernelData<F: Field> {
    kernel_nonzero_points: Vec<AffinePoint<F>>,
    a_correction_sum: F::Elem,
    b_correction_sum: F::Elem,
}

impl<F> VeluIsogeny<ShortWeierstrassCurve<F>>
where
    F: Field,
    F::Elem: Clone + Eq + Hash,
{
    fn require_non_kernel_finite_point<'a>(
        &self,
        point: &'a AffinePoint<F>,
    ) -> Result<Option<&'a AffinePoint<F>>, IsogenyError> {
        if !self.domain.contains(point) {
            return Err(IsogenyError::Curve(
                crate::elliptic_curves::CurveError::PointNotOnCurve,
            ));
        }

        if self.kernel.contains(point) {
            return Ok(None);
        }

        Ok((AffinePoint::finite_coordinates(point).is_some()).then_some(point))
    }
}

impl<F> SupportsVeluConstruction for ShortWeierstrassCurve<F>
where
    F: Field + Clone,
    F::Elem: Clone + Eq + Hash,
{
    fn velu_codomain_from_kernel(
        domain: &Self,
        kernel: &IsogenyKernel<Self>,
    ) -> Result<Self, IsogenyError> {
        VeluKernelData::from_kernel(domain, kernel).codomain_curve(domain)
    }

    fn velu_evaluate_non_kernel_point(
        isogeny: &VeluIsogeny<Self>,
        point: &Self::Point,
    ) -> Result<Self::Point, IsogenyError> {
        let finite_point = isogeny
            .require_non_kernel_finite_point(point)?
            .expect("non-kernel short-Weierstrass Vélu evaluation should be finite");
        let (x, y) = VeluKernelData::from_kernel(&isogeny.domain, &isogeny.kernel)
            .translation_correction_sums(&isogeny.domain, finite_point)?;
        isogeny.codomain.point(x, y).map_err(IsogenyError::from)
    }
}
