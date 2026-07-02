use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::traits::{
    CurveModel, CurveModelConversion, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::arb_nonsingular_general_weierstrass_curve;

use super::shared::F5;

fn reducible_transport_case() -> impl Strategy<
    Value = (
        GeneralWeierstrassCurve<F5>,
        AffinePoint<F5>,
        AffinePoint<F5>,
        u64,
    ),
> {
    arb_nonsingular_general_weierstrass_curve::<crate::fields::Fp5>(CurveStrategyConfig::default())
        .prop_flat_map(|curve| {
            let points = curve.points();
            let len = points.len();

            (
                Just(curve.clone()),
                Just(points),
                0usize..len,
                0usize..len,
                0u64..16,
            )
                .prop_map(|(curve, points, left_index, right_index, scalar)| {
                    (
                        curve,
                        points[left_index].clone(),
                        points[right_index].clone(),
                        scalar,
                    )
                })
        })
}

#[test]
fn short_reduction_preserves_classical_invariants() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert!(F5::eq(&curve.c4(), &conversion.target().c4()));
    assert!(F5::eq(&curve.c6(), &conversion.target().c6()));
    assert!(F5::eq(
        &curve.discriminant(),
        &conversion.target().discriminant()
    ));
    assert!(F5::eq(
        &curve.j_invariant(),
        &conversion.target().j_invariant()
    ));
}

#[test]
fn transport_preserves_membership_and_roundtrips_for_all_enumerated_points() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    for general_point in curve.points() {
        let short_point = conversion
            .map_source_point(&general_point)
            .expect("enumerated general point should transport to short");

        assert!(curve.contains(&general_point));
        assert!(conversion.target().contains(&short_point));
        assert_eq!(
            conversion
                .map_target_point(&short_point)
                .expect("transported short point should return to the source model"),
            general_point,
        );
    }

    for short_point in conversion.target().points() {
        let general_point = conversion
            .map_target_point(&short_point)
            .expect("enumerated short point should transport to general");

        assert!(conversion.target().contains(&short_point));
        assert!(curve.contains(&general_point));
        assert_eq!(
            conversion
                .map_source_point(&general_point)
                .expect("transported general point should return to the target model"),
            short_point,
        );
    }
}

#[test]
fn transport_preserves_addition_doubling_and_scalar_multiplication_on_a_small_curve() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");
    let scalar_limit = (curve.order() as u64) * 2 + 1;

    for left in curve.points() {
        let short_left = conversion
            .map_source_point(&left)
            .expect("enumerated point should transport to short");

        assert_eq!(
            curve.double(&left).expect("doubling should succeed"),
            conversion
                .map_target_point(
                    &conversion
                        .target()
                        .double(&short_left)
                        .expect("short doubling should succeed"),
                )
                .expect("short double should transport back"),
        );

        for scalar in 0..=scalar_limit {
            assert_eq!(
                curve
                    .mul_scalar(&left, scalar)
                    .expect("scalar multiplication should succeed"),
                conversion
                    .map_target_point(
                        &conversion
                            .target()
                            .mul_scalar(&short_left, scalar)
                            .expect("short scalar multiplication should succeed"),
                    )
                    .expect("short multiple should transport back"),
            );
        }

        for right in curve.points() {
            let short_right = conversion
                .map_source_point(&right)
                .expect("enumerated point should transport to short");

            assert_eq!(
                curve.add(&left, &right).expect("addition should succeed"),
                conversion
                    .map_target_point(
                        &conversion
                            .target()
                            .add(&short_left, &short_right)
                            .expect("short addition should succeed"),
                    )
                    .expect("short sum should transport back"),
            );
        }
    }
}

#[test]
fn transport_preserves_point_orders_and_group_order_on_a_small_curve() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(curve.order(), conversion.target().order());
    assert_eq!(curve.exponent(), conversion.target().exponent());

    for point in curve.points() {
        let short_point = conversion
            .map_source_point(&point)
            .expect("enumerated point should transport to short");

        assert_eq!(
            curve.point_order(&point),
            conversion.target().point_order(&short_point)
        );
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn property_transport_compatibility_matches_the_short_companion_on_random_reducible_curves(
        (curve, left, right, scalar) in reducible_transport_case(),
    ) {
        let conversion = curve
            .conversion_to_short_weierstrass()
            .expect("characteristic five should support the short reduction");
        let short_left = conversion
            .map_source_point(&left)
            .expect("sampled point should transport to short");
        let short_right = conversion
            .map_source_point(&right)
            .expect("sampled point should transport to short");

        prop_assert!(curve.contains(&left));
        prop_assert!(curve.contains(&right));
        prop_assert!(conversion.target().contains(&short_left));
        prop_assert!(conversion.target().contains(&short_right));
        prop_assert_eq!(curve.j_invariant(), conversion.target().j_invariant());
        prop_assert_eq!(
            curve.add(&left, &right).expect("sampled points should add"),
            conversion
                .map_target_point(
                    &conversion
                        .target()
                        .add(&short_left, &short_right)
                        .expect("short addition should succeed"),
                )
                .expect("short sum should transport back"),
        );
        prop_assert_eq!(
            curve.mul_scalar(&left, scalar)
                .expect("sampled point should multiply"),
            conversion
                .map_target_point(
                    &conversion
                        .target()
                        .mul_scalar(&short_left, scalar)
                        .expect("short scalar multiplication should succeed"),
                )
                .expect("short multiple should transport back"),
        );
    }
}
