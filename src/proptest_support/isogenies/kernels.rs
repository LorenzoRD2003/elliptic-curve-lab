use crate::elliptic_curves::{
    AffinePoint, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, ShortWeierstrassCurve,
};
use crate::fields::{Field, Fp};
use crate::isogenies::{Isogeny, VeluIsogeny};
use crate::proptest_support::combinators::same_membership_set;

pub(crate) type F41 = Fp<41>;
pub(crate) type Curve41 = ShortWeierstrassCurve<F41>;
pub(crate) type Point41 = AffinePoint<F41>;

pub(crate) fn proptest_curve41() -> Curve41 {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3))
        .expect("shared F41 curve should stay valid")
}

pub(crate) fn unique_velu_isogenies_on(curve: &Curve41) -> Vec<VeluIsogeny<Curve41>> {
    let mut kernels = Vec::<Vec<Point41>>::new();
    let mut isogenies = Vec::new();

    for point in curve.points() {
        if curve.is_identity(&point) {
            continue;
        }

        let Some(order) = curve.point_order(&point) else {
            continue;
        };
        if order <= 1 {
            continue;
        }

        let Ok(isogeny) = VeluIsogeny::from_generator(curve.clone(), point.clone()) else {
            continue;
        };
        let kernel = isogeny.kernel_points().to_vec();

        if kernels
            .iter()
            .any(|existing| same_membership_set(existing, &kernel))
        {
            continue;
        }

        kernels.push(kernel);
        isogenies.push(isogeny);
    }

    isogenies
}
