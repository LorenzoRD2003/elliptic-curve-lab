use crate::fields::traits::*;
use core::fmt;

use proptest::prelude::*;

use crate::elliptic_curves::short_weierstrass::isogenies::VeluIsogeny;
use crate::elliptic_curves::traits::{
    CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::isogenies::traits::Isogeny;
use crate::proptest_support::combinators::same_membership_set;
use crate::proptest_support::isogenies::kernels::{Curve41, Point41, proptest_curve41};

/// Small reusable Vélu kernel case over the canonical `F41` sample curve.
#[derive(Clone)]
pub struct CyclicKernelCase {
    pub curve: Curve41,
    pub generator: Point41,
    pub kernel_point: Point41,
    pub sample_point: Point41,
    pub coset_point: Point41,
    pub isogeny: VeluIsogeny<Curve41>,
}

impl fmt::Debug for CyclicKernelCase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CyclicKernelCase")
            .field("generator", &self.generator)
            .field("kernel_point", &self.kernel_point)
            .field("sample_point", &self.sample_point)
            .field("coset_point", &self.coset_point)
            .field("kernel_size", &self.isogeny.kernel_points().len())
            .finish()
    }
}

/// Returns one tiny cyclic-kernel Vélu fixture over the canonical `F41` curve.
pub fn arb_cyclic_kernel_case() -> BoxedStrategy<CyclicKernelCase> {
    prop::sample::select(all_cyclic_kernel_cases()).boxed()
}

pub(crate) fn all_cyclic_kernel_cases() -> Vec<CyclicKernelCase> {
    let curve = proptest_curve41();
    let all_points = curve.points();
    let mut kernels = Vec::<Vec<Point41>>::new();
    let mut cases = Vec::new();

    for generator in &all_points {
        if curve.is_identity(generator) {
            continue;
        }

        let Some(order) = curve.point_order(generator) else {
            continue;
        };
        if order <= 1 {
            continue;
        }

        let Ok(isogeny) = VeluIsogeny::from_generator(curve.clone(), generator.clone()) else {
            continue;
        };
        let kernel_points = isogeny.kernel_points().to_vec();
        if kernel_points.len() == all_points.len() {
            continue;
        }
        if kernels
            .iter()
            .any(|existing| same_membership_set(existing, &kernel_points))
        {
            continue;
        }
        kernels.push(kernel_points.clone());

        let kernel_point = kernel_points
            .iter()
            .find(|point| !curve.is_identity(point))
            .expect("non-trivial cyclic kernel should contain a non-identity point")
            .clone();
        let sample_point = all_points
            .iter()
            .find(|point| !kernel_points.contains(point))
            .expect("sample curve should have points outside the kernel")
            .clone();
        let coset_point = curve
            .add(&sample_point, &kernel_point)
            .expect("kernel translation should stay on the curve");

        cases.push(CyclicKernelCase {
            curve: curve.clone(),
            generator: generator.clone(),
            kernel_point,
            sample_point,
            coset_point,
            isogeny,
        });
    }

    cases
}

pub(crate) fn touch_cyclic_kernel_case_fields() {
    let _ = |case: CyclicKernelCase| {
        let _ = case.curve;
        let _ = case.generator;
        let _ = case.kernel_point;
        let _ = case.sample_point;
        let _ = case.coset_point;
        let _ = case.isogeny;
    };
}
