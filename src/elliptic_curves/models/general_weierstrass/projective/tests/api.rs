use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::traits::{
    AffineCurveModel, CurveModelConversion, EnumerableCurveModel, GroupCurveModel,
    HasProjectiveModel, ProjectiveGroupCurveModel,
};
use crate::elliptic_curves::{AffinePoint, GeneralWeierstrassCurve, ProjectivePoint};
use crate::proptest_support::{
    config::CurveStrategyConfig,
    elliptic_curves::{
        arb_general_weierstrass_projective_equivalence_class,
        arb_general_weierstrass_projective_pair, arb_general_weierstrass_projective_point,
        rescale_projective_point,
    },
};

use crate::elliptic_curves::models::general_weierstrass::projective::{
    GeneralWeierstrassProjectiveOperationCost, GeneralWeierstrassProjectiveOperationKind,
};

type F2 = crate::fields::Fp2;
type F3 = crate::fields::Fp3;
type F5 = crate::fields::Fp5;

fn curve_f5() -> GeneralWeierstrassCurve<F5> {
    GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
        .expect("non-singular curve")
}

fn scaled_projective(point: &AffinePoint<F5>, scale: i64) -> ProjectivePoint<F5> {
    match point {
        AffinePoint::Infinity => ProjectivePoint::Infinity,
        AffinePoint::Finite { x, y } => {
            let lambda = F5::from_i64(scale);
            ProjectivePoint::new(F5::mul(x, &lambda), F5::mul(y, &lambda), lambda)
        }
    }
}

#[test]
fn general_projective_roundtrip_recovers_the_same_affine_point() {
    let curve = curve_f5();
    let point = curve
        .point(F5::from_i64(2), F5::one())
        .expect("point should lie on the curve");
    let projective = scaled_projective(&point, 2);

    assert!(curve.is_projective_point_on_curve(&projective));
    assert_eq!(curve.to_affine_projective(&projective), Ok(point));
    assert_eq!(
        projective.normalize(),
        curve.to_projective(
            &curve
                .point(F5::from_i64(2), F5::one())
                .expect("point should lie on the curve")
        )
    );
}

#[test]
fn general_projective_negation_uses_the_model_specific_involution() {
    let curve = curve_f5();
    let point = curve
        .point(F5::zero(), F5::from_i64(4))
        .expect("point should lie on the curve");
    let projective = scaled_projective(&point, 3);

    assert_eq!(
        curve.to_affine_projective(&curve.neg_projective(&projective)),
        Ok(curve.neg(&point))
    );
    assert_ne!(curve.neg(&point), point.neg());
}

#[test]
fn general_projective_group_operations_match_the_affine_route_on_small_examples() {
    let curve = curve_f5();
    let left = curve
        .point(F5::zero(), F5::zero())
        .expect("left point should lie on the curve");
    let right = curve
        .point(F5::from_i64(2), F5::one())
        .expect("right point should lie on the curve");
    let left_projective = scaled_projective(&left, 2);
    let right_projective = scaled_projective(&right, 4);

    assert_eq!(
        curve.to_affine_projective(
            &curve
                .add_projective(&left_projective, &right_projective)
                .expect("projective sum should succeed")
        ),
        Ok(curve.add(&left, &right).expect("affine sum should succeed"))
    );
    assert_eq!(
        curve.to_affine_projective(
            &curve
                .double_projective(&left_projective)
                .expect("projective doubling should succeed")
        ),
        Ok(curve.double(&left).expect("affine double should succeed"))
    );
    assert_eq!(
        curve.to_affine_projective(
            &curve
                .mixed_add_projective(&left_projective, &right)
                .expect("projective mixed addition should succeed")
        ),
        Ok(curve.add(&left, &right).expect("affine sum should succeed"))
    );
    assert_eq!(
        curve.to_affine_projective(
            &curve
                .mul_scalar_projective(&left_projective, 2)
                .expect("projective scalar multiplication should succeed")
        ),
        Ok(curve
            .mul_scalar(&left, 2)
            .expect("affine multiple should succeed"))
    );
}

#[test]
fn general_projective_operations_respect_characteristics_two_and_three() {
    let curve_two =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let point_two = curve_two
        .point(F2::one(), F2::zero())
        .expect("point should lie on the curve");
    let point_two_projective = curve_two
        .to_projective(&point_two)
        .expect("point should lift");

    assert_eq!(
        curve_two.to_affine_projective(
            &curve_two
                .double_projective(&point_two_projective)
                .expect("projective doubling should succeed")
        ),
        Ok(curve_two
            .double(&point_two)
            .expect("affine double should succeed"))
    );

    let curve_three =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let point_three = curve_three
        .point(F3::zero(), F3::zero())
        .expect("point should lie on the curve");
    let point_three_projective = curve_three
        .to_projective(&point_three)
        .expect("point should lift");

    assert_eq!(
        curve_three.to_affine_projective(
            &curve_three
                .mul_scalar_projective(&point_three_projective, 2)
                .expect("projective scalar multiplication should succeed")
        ),
        Ok(curve_three
            .mul_scalar(&point_three, 2)
            .expect("affine multiple should succeed"))
    );
}

#[test]
fn general_projective_cost_reports_remain_stable_for_the_native_projective_story() {
    assert_eq!(
        GeneralWeierstrassProjectiveOperationCost::for_kind(
            GeneralWeierstrassProjectiveOperationKind::Neg,
        )
        .representation_cost()
        .multiplications(),
        2
    );
    assert_eq!(
        GeneralWeierstrassProjectiveOperationCost::for_kind(
            GeneralWeierstrassProjectiveOperationKind::Add,
        )
        .affine_additions(),
        0
    );
    assert_eq!(
        GeneralWeierstrassProjectiveOperationCost::for_scalar_mul(13)
            .representation_cost()
            .inversions(),
        0
    );
}

#[test]
fn characteristic_three_scaled_representatives_roundtrip_exhaustively() {
    let scales = [F3::one(), F3::from_i64(2)];

    for a1 in 0..3 {
        for a2 in 0..3 {
            for a3 in 0..3 {
                for a4 in 0..3 {
                    for a6 in 0..3 {
                        let Ok(curve) = GeneralWeierstrassCurve::<F3>::new(
                            F3::from_i64(a1),
                            F3::from_i64(a2),
                            F3::from_i64(a3),
                            F3::from_i64(a4),
                            F3::from_i64(a6),
                        ) else {
                            continue;
                        };

                        for point in curve.points() {
                            for scale in scales {
                                let projective = match &point {
                                    AffinePoint::Infinity => ProjectivePoint::Infinity,
                                    AffinePoint::Finite { x, y } => ProjectivePoint::new(
                                        F3::mul(x, &scale),
                                        F3::mul(y, &scale),
                                        scale,
                                    ),
                                };
                                assert!(curve.is_projective_point_on_curve(&projective));
                                assert_eq!(
                                    curve.to_affine_projective(&projective),
                                    Ok(point.clone())
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_general_projective_roundtrip_recovers_the_same_affine_point(
        (curve, affine, projective) in arb_general_weierstrass_projective_point::<crate::fields::Fp5>(CurveStrategyConfig::default()),
    ) {
        prop_assert!(curve.is_projective_point_on_curve(&projective));
        prop_assert_eq!(curve.to_affine_projective(&projective), Ok(affine.clone()));
        prop_assert_eq!(projective.normalize(), curve.to_projective(&affine));
    }

    #[test]
    fn property_equivalent_general_projective_representatives_share_affine_image_and_normalized_chart_in_characteristic_five(
        (curve, affine, left_projective, right_projective) in
            arb_general_weierstrass_projective_equivalence_class::<crate::fields::Fp5>(CurveStrategyConfig::default()),
    ) {
        prop_assert!(curve.is_projective_point_on_curve(&left_projective));
        prop_assert!(curve.is_projective_point_on_curve(&right_projective));
        prop_assert_eq!(curve.to_affine_projective(&left_projective), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&right_projective), Ok(affine.clone()));
        prop_assert_eq!(left_projective.normalize(), right_projective.normalize());
        prop_assert_eq!(left_projective.normalize(), curve.to_projective(&affine));
    }

    #[test]
    fn property_general_projective_group_operations_match_the_affine_route_in_characteristic_two(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_general_weierstrass_projective_pair::<crate::fields::Fp2>(CurveStrategyConfig::default()),
        scalar in 0u64..16,
    ) {
        let projective_sum = curve.add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let projective_double = curve.double_projective(&left_projective)
            .expect("projective doubling should succeed");
        let projective_multiple = curve.mul_scalar_projective(&left_projective, scalar)
            .expect("projective scalar multiplication should succeed");

        prop_assert_eq!(curve.to_affine_projective(&projective_sum), Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_double), Ok(curve.double(&left_affine).expect("affine double should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_multiple), Ok(curve.mul_scalar(&left_affine, scalar).expect("affine multiple should succeed")));
    }

    #[test]
    fn property_general_projective_identity_inverse_and_doubling_laws_in_characteristic_two(
        (curve, affine, projective) in arb_general_weierstrass_projective_point::<crate::fields::Fp2>(CurveStrategyConfig::default()),
    ) {
        let identity = curve.projective_identity();
        let negated = curve.neg_projective(&projective);
        let left_identity_sum = curve.add_projective(&projective, &identity)
            .expect("adding the projective identity should succeed");
        let right_identity_sum = curve.add_projective(&identity, &projective)
            .expect("adding the projective identity should succeed");
        let inverse_sum = curve.add_projective(&projective, &negated)
            .expect("adding a projective inverse should succeed");
        let doubled = curve.double_projective(&projective)
            .expect("projective doubling should succeed");
        let self_sum = curve.add_projective(&projective, &projective)
            .expect("projective self-addition should succeed");

        prop_assert_eq!(curve.to_affine_projective(&left_identity_sum), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&right_identity_sum), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&inverse_sum), Ok(AffinePoint::Infinity));
        prop_assert_eq!(curve.to_affine_projective(&doubled), curve.to_affine_projective(&self_sum));
    }

    #[test]
    fn property_general_projective_group_operations_match_the_affine_route_in_characteristic_three(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_general_weierstrass_projective_pair::<crate::fields::Fp3>(CurveStrategyConfig::default()),
        scalar in 0u64..16,
    ) {
        let projective_sum = curve.add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let projective_double = curve.double_projective(&left_projective)
            .expect("projective doubling should succeed");
        let projective_multiple = curve.mul_scalar_projective(&left_projective, scalar)
            .expect("projective scalar multiplication should succeed");

        prop_assert_eq!(curve.to_affine_projective(&projective_sum), Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_double), Ok(curve.double(&left_affine).expect("affine double should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_multiple), Ok(curve.mul_scalar(&left_affine, scalar).expect("affine multiple should succeed")));
    }

    #[test]
    fn property_general_projective_identity_inverse_and_doubling_laws_in_characteristic_three(
        (curve, affine, projective) in arb_general_weierstrass_projective_point::<crate::fields::Fp3>(CurveStrategyConfig::default()),
    ) {
        let identity = curve.projective_identity();
        let negated = curve.neg_projective(&projective);
        let left_identity_sum = curve.add_projective(&projective, &identity)
            .expect("adding the projective identity should succeed");
        let right_identity_sum = curve.add_projective(&identity, &projective)
            .expect("adding the projective identity should succeed");
        let inverse_sum = curve.add_projective(&projective, &negated)
            .expect("adding a projective inverse should succeed");
        let doubled = curve.double_projective(&projective)
            .expect("projective doubling should succeed");
        let self_sum = curve.add_projective(&projective, &projective)
            .expect("projective self-addition should succeed");

        prop_assert_eq!(curve.to_affine_projective(&left_identity_sum), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&right_identity_sum), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&inverse_sum), Ok(AffinePoint::Infinity));
        prop_assert_eq!(curve.to_affine_projective(&doubled), curve.to_affine_projective(&self_sum));
    }

    #[test]
    fn property_general_projective_group_operations_match_the_affine_route_and_short_companion_in_characteristic_five(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_general_weierstrass_projective_pair::<crate::fields::Fp5>(CurveStrategyConfig::default()),
        scalar in 0u64..16,
    ) {
        let conversion = curve
            .conversion_to_short_weierstrass()
            .expect("characteristic five should support the short reduction");
        let short_left = conversion.map_source_point(&left_affine)
            .expect("left point should transport to short");
        let short_right = conversion.map_source_point(&right_affine)
            .expect("right point should transport to short");

        let projective_sum = curve.add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let projective_double = curve.double_projective(&left_projective)
            .expect("projective doubling should succeed");
        let projective_mixed = curve.mixed_add_projective(&left_projective, &right_affine)
            .expect("projective mixed addition should succeed");
        let projective_multiple = curve.mul_scalar_projective(&left_projective, scalar)
            .expect("projective scalar multiplication should succeed");

        prop_assert_eq!(curve.to_affine_projective(&projective_sum), Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_double), Ok(curve.double(&left_affine).expect("affine double should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_mixed), Ok(curve.add(&left_affine, &right_affine).expect("affine mixed sum should succeed")));
        prop_assert_eq!(curve.to_affine_projective(&projective_multiple), Ok(curve.mul_scalar(&left_affine, scalar).expect("affine multiple should succeed")));

        prop_assert_eq!(
            curve.to_affine_projective(&projective_sum),
            Ok(conversion.map_target_point(&conversion.target().add(&short_left, &short_right).expect("short sum should succeed")).expect("short sum should transport back"))
        );
        prop_assert_eq!(
            curve.to_affine_projective(&projective_double),
            Ok(conversion.map_target_point(&conversion.target().double(&short_left).expect("short double should succeed")).expect("short double should transport back"))
        );
        prop_assert_eq!(
            curve.to_affine_projective(&projective_multiple),
            Ok(conversion.map_target_point(&conversion.target().mul_scalar(&short_left, scalar).expect("short multiple should succeed")).expect("short multiple should transport back"))
        );
    }

    #[test]
    fn property_general_projective_identity_inverse_doubling_rescaling_and_normalization_in_characteristic_five(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_general_weierstrass_projective_pair::<crate::fields::Fp5>(CurveStrategyConfig::default()),
        left_scale in 1u64..5,
        right_scale in 1u64..5,
        scalar in 0u64..16,
    ) {
        let identity = curve.projective_identity();
        let negated = curve.neg_projective(&left_projective);
        let left_identity_sum = curve.add_projective(&left_projective, &identity)
            .expect("adding the projective identity should succeed");
        let right_identity_sum = curve.add_projective(&identity, &left_projective)
            .expect("adding the projective identity should succeed");
        let inverse_sum = curve.add_projective(&left_projective, &negated)
            .expect("adding a projective inverse should succeed");
        let doubled = curve.double_projective(&left_projective)
            .expect("projective doubling should succeed");
        let self_sum = curve.add_projective(&left_projective, &left_projective)
            .expect("projective self-addition should succeed");

        let left_rescaled = rescale_projective_point::<crate::fields::Fp5>(&left_projective, &crate::fields::Fp5::from_i64(left_scale));
        let right_rescaled = rescale_projective_point::<crate::fields::Fp5>(&right_projective, &crate::fields::Fp5::from_i64(right_scale));

        let sum = curve.add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let rescaled_sum = curve.add_projective(&left_rescaled, &right_rescaled)
            .expect("projective addition should succeed");
        let normalized_sum = curve.add_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
            &right_projective.normalize().expect("valid projective point should normalize"),
        ).expect("projective addition should succeed");

        let rescaled_doubled = curve.double_projective(&left_rescaled)
            .expect("projective doubling should succeed");
        let normalized_doubled = curve.double_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
        ).expect("projective doubling should succeed");

        let multiple = curve.mul_scalar_projective(&left_projective, scalar)
            .expect("projective scalar multiplication should succeed");
        let rescaled_multiple = curve.mul_scalar_projective(&left_rescaled, scalar)
            .expect("projective scalar multiplication should succeed");
        let normalized_multiple = curve.mul_scalar_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
            scalar,
        ).expect("projective scalar multiplication should succeed");

        let mixed_sum = curve.mixed_add_projective(&left_projective, &right_affine)
            .expect("projective mixed addition should succeed");
        let normalized_mixed_sum = curve.mixed_add_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
            &right_affine,
        ).expect("projective mixed addition should succeed");

        prop_assert_eq!(curve.to_affine_projective(&left_identity_sum), Ok(left_affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&right_identity_sum), Ok(left_affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&inverse_sum), Ok(AffinePoint::Infinity));
        prop_assert_eq!(curve.to_affine_projective(&doubled), curve.to_affine_projective(&self_sum));

        prop_assert_eq!(curve.to_affine_projective(&sum), curve.to_affine_projective(&rescaled_sum));
        prop_assert_eq!(curve.to_affine_projective(&sum), curve.to_affine_projective(&normalized_sum));
        prop_assert_eq!(curve.to_affine_projective(&doubled), curve.to_affine_projective(&rescaled_doubled));
        prop_assert_eq!(curve.to_affine_projective(&doubled), curve.to_affine_projective(&normalized_doubled));
        prop_assert_eq!(curve.to_affine_projective(&multiple), curve.to_affine_projective(&rescaled_multiple));
        prop_assert_eq!(curve.to_affine_projective(&multiple), curve.to_affine_projective(&normalized_multiple));
        prop_assert_eq!(curve.to_affine_projective(&mixed_sum), curve.to_affine_projective(&normalized_mixed_sum));
        prop_assert_eq!(
            curve.to_affine_projective(&sum),
            Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed"))
        );
    }
}
