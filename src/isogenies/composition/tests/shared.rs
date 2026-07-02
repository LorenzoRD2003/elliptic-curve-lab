use crate::elliptic_curves::{
    AffinePoint, ShortWeierstrassCurve,
    short_weierstrass::isomorphisms::ShortWeierstrassIsomorphism,
    traits::{
        AffineCurveModel, CurveIsomorphism, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel,
    },
};
use crate::isogenies::{
    error::IsogenyError,
    kernel::{KernelDescription, ReducedKernelDescription},
    traits::Isogeny,
    velu::VeluIsogeny,
};

pub(super) type F41 = crate::fields::Fp41;
pub(super) type Curve = ShortWeierstrassCurve<F41>;
pub(super) type Velu = VeluIsogeny<Curve>;

pub(super) fn curve_a() -> Curve {
    Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
}

pub(super) fn first_generator(curve: &Curve) -> AffinePoint<F41> {
    curve
        .point(F41::from_i64(40), F41::from_i64(0))
        .expect("known two-torsion point should lie on the curve")
}

pub(super) fn small_nontrivial_generator(curve: &Curve) -> AffinePoint<F41> {
    curve
        .point_orders()
        .into_iter()
        .filter(|(point, order)| *order > 1 && !curve.is_identity(point))
        .min_by_key(|(_, order)| *order)
        .map(|(point, _)| point)
        .expect("small sample curve should have a non-trivial point")
}

pub(super) fn first_isogeny() -> Velu {
    let domain = curve_a();
    VeluIsogeny::from_generator(domain.clone(), first_generator(&domain))
        .expect("first sample Vélu isogeny should build")
}

pub(super) fn second_isogeny(domain: &Curve) -> Velu {
    VeluIsogeny::from_generator(domain.clone(), small_nontrivial_generator(domain))
        .expect("second sample Vélu isogeny should build")
}

pub(super) fn bridged_second_isogeny(first: &Velu) -> (ShortWeierstrassIsomorphism<F41>, Velu) {
    let bridge = ShortWeierstrassIsomorphism::new(first.codomain().clone(), F41::from_i64(3))
        .expect("sample bridge isomorphism should build");
    let generator = small_nontrivial_generator(first.codomain());
    let transported_generator = bridge
        .evaluate(&generator)
        .expect("bridge should transport the sample generator");
    let second = VeluIsogeny::from_generator(bridge.codomain().clone(), transported_generator)
        .expect("bridged second sample Vélu isogeny should build");

    (bridge, second)
}

pub(super) fn first_non_kernel_point(isogeny: &Velu) -> AffinePoint<F41> {
    isogeny
        .domain()
        .points()
        .into_iter()
        .find(|point| !isogeny.kernel_points().contains(point))
        .expect("sample Vélu isogeny should have at least one point outside its kernel")
}

pub(super) struct BrokenMiddleImageIsogeny {
    pub(super) inner: Velu,
    pub(super) broken_point: AffinePoint<F41>,
}

impl Isogeny<Curve, Curve> for BrokenMiddleImageIsogeny {
    fn domain(&self) -> &Curve {
        self.inner.domain()
    }

    fn codomain(&self) -> &Curve {
        self.inner.codomain()
    }

    fn degree(&self) -> usize {
        self.inner.degree()
    }

    fn evaluate(
        &self,
        point: &<Curve as CurveModel>::Point,
    ) -> Result<<Curve as CurveModel>::Point, IsogenyError> {
        if *point == self.broken_point {
            return Ok(AffinePoint::new(F41::from_i64(2), F41::from_i64(2)));
        }

        self.inner.evaluate(point)
    }

    fn kernel_description(&self) -> KernelDescription<Curve> {
        KernelDescription::Reduced(
            ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints {
                points: self.inner.kernel_points(),
                degree: self.inner.degree(),
            },
        )
    }
}
