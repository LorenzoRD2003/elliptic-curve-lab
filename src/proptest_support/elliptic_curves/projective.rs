use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::general_weierstrass::GeneralWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};
use crate::elliptic_curves::{AffinePoint, ProjectivePoint};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::{
    arb_curve_and_point, arb_general_weierstrass_curve_and_point,
    general_weierstrass::arb_nonsingular_general_weierstrass_curve,
    short_weierstrass::arb_nonsingular_curve,
};

fn scaled_projective_point<F>(
    point: &AffinePoint<F>,
    scale: &<F as Field>::Elem,
) -> ProjectivePoint<F>
where
    F: Field,
{
    match point {
        AffinePoint::Infinity => ProjectivePoint::Infinity,
        AffinePoint::Finite { x, y } => {
            ProjectivePoint::new(F::mul(x, scale), F::mul(y, scale), scale.clone())
        }
    }
}

/// Returns one projective representative in the same homogeneous equivalence
/// class by multiplying every finite coordinate by the same nonzero field
/// element.
pub fn rescale_projective_point<F>(
    point: &ProjectivePoint<F>,
    scale: &<F as Field>::Elem,
) -> ProjectivePoint<F>
where
    F: Field,
{
    match point {
        ProjectivePoint::Infinity => ProjectivePoint::Infinity,
        ProjectivePoint::Finite { x, y, z } => {
            ProjectivePoint::new(F::mul(x, scale), F::mul(y, scale), F::mul(z, scale))
        }
    }
}

/// Returns a short-Weierstrass curve together with one affine/projective pair
/// representing the same rational point.
pub fn arb_short_weierstrass_projective_point<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(ShortWeierstrassCurve<F>, AffinePoint<F>, ProjectivePoint<F>)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_curve_and_point::<F>(config)
        .prop_flat_map(move |(curve, point)| {
            let point_for_projection = point.clone();
            (
                Just(curve),
                Just(point),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(move |(curve, point, scale)| {
                    (
                        curve,
                        point,
                        scaled_projective_point(&point_for_projection, &scale),
                    )
                })
        })
        .boxed()
}

/// Returns a short-Weierstrass curve together with one affine point and two
/// projective representatives from the same homogeneous equivalence class.
pub fn arb_short_weierstrass_projective_equivalence_class<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(
    ShortWeierstrassCurve<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
    ProjectivePoint<F>,
)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_curve_and_point::<F>(config)
        .prop_flat_map(move |(curve, point)| {
            let point_for_left = point.clone();
            let point_for_right = point.clone();
            (
                Just(curve),
                Just(point),
                prop::sample::select(nonzero_scales.clone()),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(move |(curve, point, left_scale, right_scale)| {
                    (
                        curve,
                        point,
                        scaled_projective_point(&point_for_left, &left_scale),
                        scaled_projective_point(&point_for_right, &right_scale),
                    )
                })
        })
        .boxed()
}

/// Returns one short-Weierstrass curve together with two affine/projective
/// pairs on that same curve.
pub fn arb_short_weierstrass_projective_pair<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(
    ShortWeierstrassCurve<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_nonsingular_curve::<F>(config)
        .prop_flat_map(move |curve| {
            let points = curve.points();
            let point_count = points.len();

            (
                Just(curve.clone()),
                Just(points),
                0usize..point_count,
                0usize..point_count,
                prop::sample::select(nonzero_scales.clone()),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(
                    |(curve, points, left_index, right_index, left_scale, right_scale)| {
                        let left = points[left_index].clone();
                        let right = points[right_index].clone();
                        let left_projective = scaled_projective_point(&left, &left_scale);
                        let right_projective = scaled_projective_point(&right, &right_scale);

                        (curve, left, left_projective, right, right_projective)
                    },
                )
        })
        .boxed()
}

/// Returns a general-Weierstrass curve together with one affine/projective pair
/// representing the same rational point.
pub fn arb_general_weierstrass_projective_point<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(
    GeneralWeierstrassCurve<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
)>
where
    F: EnumerableFiniteField + 'static,
    F::Elem: 'static,
    GeneralWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_general_weierstrass_curve_and_point::<F>(config)
        .prop_flat_map(move |(curve, point)| {
            let point_for_projection = point.clone();
            (
                Just(curve),
                Just(point),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(move |(curve, point, scale)| {
                    (
                        curve,
                        point,
                        scaled_projective_point(&point_for_projection, &scale),
                    )
                })
        })
        .boxed()
}

/// Returns a general-Weierstrass curve together with one affine point and two
/// projective representatives from the same homogeneous equivalence class.
pub fn arb_general_weierstrass_projective_equivalence_class<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(
    GeneralWeierstrassCurve<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
    ProjectivePoint<F>,
)>
where
    F: EnumerableFiniteField + 'static,
    F::Elem: 'static,
    GeneralWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_general_weierstrass_curve_and_point::<F>(config)
        .prop_flat_map(move |(curve, point)| {
            let point_for_left = point.clone();
            let point_for_right = point.clone();
            (
                Just(curve),
                Just(point),
                prop::sample::select(nonzero_scales.clone()),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(move |(curve, point, left_scale, right_scale)| {
                    (
                        curve,
                        point,
                        scaled_projective_point(&point_for_left, &left_scale),
                        scaled_projective_point(&point_for_right, &right_scale),
                    )
                })
        })
        .boxed()
}

/// Returns one general-Weierstrass curve together with two affine/projective
/// pairs on that same curve.
pub fn arb_general_weierstrass_projective_pair<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<(
    GeneralWeierstrassCurve<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
    AffinePoint<F>,
    ProjectivePoint<F>,
)>
where
    F: EnumerableFiniteField + 'static,
    F::Elem: 'static,
    GeneralWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    let nonzero_scales = F::elements()
        .into_iter()
        .filter(|value| !F::is_zero(value))
        .collect::<Vec<_>>();

    arb_nonsingular_general_weierstrass_curve::<F>(config)
        .prop_flat_map(move |curve| {
            let points = curve.points();
            let point_count = points.len();

            (
                Just(curve.clone()),
                Just(points),
                0usize..point_count,
                0usize..point_count,
                prop::sample::select(nonzero_scales.clone()),
                prop::sample::select(nonzero_scales.clone()),
            )
                .prop_map(
                    |(curve, points, left_index, right_index, left_scale, right_scale)| {
                        let left = points[left_index].clone();
                        let right = points[right_index].clone();
                        let left_projective = scaled_projective_point(&left, &left_scale);
                        let right_projective = scaled_projective_point(&right, &right_scale);

                        (curve, left, left_projective, right, right_projective)
                    },
                )
        })
        .boxed()
}
