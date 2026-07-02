use crate::fields::traits::*;
use proptest::prelude::*;

use crate::elliptic_curves::{
    AffinePoint, ProjectivePoint, ShortWeierstrassCurve,
    short_weierstrass::projective::{
        ShortWeierstrassProjectiveOperationCost, ShortWeierstrassProjectiveOperationKind,
    },
    traits::{
        AffineCurveModel, EnumerableCurveModel, GroupCurveModel, HasProjectiveModel,
        ProjectiveGroupCurveModel,
    },
};
use crate::proptest_support::{
    config::CurveStrategyConfig,
    elliptic_curves::{
        arb_short_weierstrass_projective_equivalence_class, arb_short_weierstrass_projective_pair,
        arb_short_weierstrass_projective_point, rescale_projective_point,
    },
};

type F7 = crate::fields::Fp7;

fn f7_curve() -> ShortWeierstrassCurve<F7> {
    ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
}

fn scaled_projective(point: &AffinePoint<F7>, scale: i64) -> ProjectivePoint<F7> {
    match point {
        AffinePoint::Infinity => ProjectivePoint::Infinity,
        AffinePoint::Finite { x, y } => {
            let lambda = F7::from_i64(scale);
            ProjectivePoint::new(F7::mul(x, &lambda), F7::mul(y, &lambda), lambda)
        }
    }
}

#[test]
fn projective_membership_accepts_scaled_representatives_of_curve_points() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let projective = scaled_projective(&point, 3);

    assert!(curve.is_projective_point_on_curve(&projective));
    assert_eq!(curve.to_affine_projective(&projective), Ok(point));
}

#[test]
fn projective_negation_preserves_the_same_affine_inverse() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let projective = scaled_projective(&point, 3);

    assert_eq!(
        curve.to_affine_projective(&curve.neg_projective(&projective)),
        Ok(curve.neg(&point))
    );
}

#[test]
fn projective_group_operations_match_the_affine_group_law_on_known_examples() {
    let curve = f7_curve();
    let left = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let right = curve
        .point(F7::from_i64(3), F7::from_i64(1))
        .expect("point should lie on the curve");
    let left_projective = scaled_projective(&left, 3);
    let right_projective = scaled_projective(&right, 5);

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
                .expect("mixed addition should succeed")
        ),
        Ok(curve.add(&left, &right).expect("affine sum should succeed"))
    );
}

#[test]
fn projective_scalar_multiplication_matches_the_affine_route() {
    let curve = f7_curve();
    let point = curve
        .point(F7::from_i64(2), F7::from_i64(1))
        .expect("point should lie on the curve");
    let projective = scaled_projective(&point, 3);

    assert_eq!(
        curve.to_affine_projective(
            &curve
                .mul_scalar_projective(&projective, 6)
                .expect("projective scalar multiplication should succeed")
        ),
        Ok(curve
            .mul_scalar(&point, 6)
            .expect("affine multiple should succeed"))
    );
}

#[test]
fn projective_normalization_rejects_scaled_triples_with_zero_z() {
    let curve = f7_curve();
    let invalid = ProjectivePoint::<F7>::new(F7::from_i64(2), F7::from_i64(1), F7::zero());

    assert!(!curve.is_projective_point_on_curve(&invalid));
    assert!(curve.to_affine_projective(&invalid).is_err());
    assert!(invalid.normalize().is_err());
}

#[test]
fn educational_cost_reports_remain_stable_for_the_baseline_bridge_story() {
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::FromAffine,
        )
        .representation_cost()
        .inversions(),
        0
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::ToAffine,
        )
        .representation_cost()
        .inversions(),
        1
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Normalize,
        )
        .representation_cost()
        .multiplications(),
        2
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Neg,
        )
        .affine_additions(),
        0
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Add,
        )
        .affine_additions(),
        0
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::Double,
        )
        .affine_doublings(),
        0
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_kind(
            ShortWeierstrassProjectiveOperationKind::MixedAdd,
        )
        .representation_cost()
        .inversions(),
        0
    );
    assert_eq!(
        ShortWeierstrassProjectiveOperationCost::for_scalar_mul(13)
            .representation_cost()
            .inversions(),
        0
    );
}

#[test]
fn exhaustive_scaled_representatives_roundtrip_over_a_small_curve() {
    let curve = f7_curve();
    let scales = [1_i64, 2, 3, 4, 5, 6];

    for point in curve.points() {
        for scale in scales {
            let projective = scaled_projective(&point, scale);
            assert!(curve.is_projective_point_on_curve(&projective));
            assert_eq!(curve.to_affine_projective(&projective), Ok(point.clone()));
            assert_eq!(projective.normalize(), curve.to_projective(&point));
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_projective_roundtrip_recovers_the_same_affine_point(
        (curve, affine, projective) in arb_short_weierstrass_projective_point::<crate::fields::Fp17>(CurveStrategyConfig::default()),
    ) {
        prop_assert!(curve.is_projective_point_on_curve(&projective));
        prop_assert_eq!(curve.to_affine_projective(&projective), Ok(affine.clone()));
        prop_assert_eq!(projective.normalize(), curve.to_projective(&affine));
    }

    #[test]
    fn property_equivalent_short_projective_representatives_share_affine_image_and_normalized_chart(
        (curve, affine, left_projective, right_projective) in
            arb_short_weierstrass_projective_equivalence_class::<crate::fields::Fp17>(CurveStrategyConfig::default()),
    ) {
        prop_assert!(curve.is_projective_point_on_curve(&left_projective));
        prop_assert!(curve.is_projective_point_on_curve(&right_projective));
        prop_assert_eq!(curve.to_affine_projective(&left_projective), Ok(affine.clone()));
        prop_assert_eq!(curve.to_affine_projective(&right_projective), Ok(affine.clone()));
        prop_assert_eq!(left_projective.normalize(), right_projective.normalize());
        prop_assert_eq!(left_projective.normalize(), curve.to_projective(&affine));
    }

    #[test]
    fn property_projective_group_operations_match_the_affine_route(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_short_weierstrass_projective_pair::<crate::fields::Fp17>(CurveStrategyConfig::default()),
        scalar in 0u64..16,
    ) {
        let projective_sum = curve
            .add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let projective_double = curve
            .double_projective(&left_projective)
            .expect("projective doubling should succeed");
        let projective_mixed = curve
            .mixed_add_projective(&left_projective, &right_affine)
            .expect("mixed addition should succeed");
        let projective_multiple = curve
            .mul_scalar_projective(&left_projective, scalar)
            .expect("projective scalar multiplication should succeed");

        prop_assert_eq!(
            curve.to_affine_projective(&projective_sum),
            Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed"))
        );
        prop_assert_eq!(
            curve.to_affine_projective(&projective_double),
            Ok(curve.double(&left_affine).expect("affine double should succeed"))
        );
        prop_assert_eq!(
            curve.to_affine_projective(&projective_mixed),
            Ok(curve.add(&left_affine, &right_affine).expect("affine sum should succeed"))
        );
        prop_assert_eq!(
            curve.to_affine_projective(&projective_multiple),
            Ok(curve.mul_scalar(&left_affine, scalar).expect("affine scalar multiplication should succeed"))
        );
    }

    #[test]
    fn property_projective_identity_inverse_and_doubling_laws(
        (curve, affine, projective) in arb_short_weierstrass_projective_point::<crate::fields::Fp17>(CurveStrategyConfig::default()),
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
        prop_assert_eq!(
            curve.to_affine_projective(&doubled),
            curve.to_affine_projective(&self_sum),
        );
    }

    #[test]
    fn property_projective_operations_are_invariant_under_rescaling_and_normalization(
        (curve, left_affine, left_projective, right_affine, right_projective) in
            arb_short_weierstrass_projective_pair::<crate::fields::Fp17>(CurveStrategyConfig::default()),
        left_scale in 1u64..17,
        right_scale in 1u64..17,
        scalar in 0u64..16,
    ) {
        let left_rescaled = rescale_projective_point::<crate::fields::Fp17>(&left_projective, &crate::fields::Fp17::from_i64(left_scale));
        let right_rescaled = rescale_projective_point::<crate::fields::Fp17>(&right_projective, &crate::fields::Fp17::from_i64(right_scale));

        let sum = curve.add_projective(&left_projective, &right_projective)
            .expect("projective addition should succeed");
        let rescaled_sum = curve.add_projective(&left_rescaled, &right_rescaled)
            .expect("projective addition should succeed");
        let normalized_sum = curve.add_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
            &right_projective.normalize().expect("valid projective point should normalize"),
        ).expect("projective addition should succeed");

        let doubled = curve.double_projective(&left_projective)
            .expect("projective doubling should succeed");
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
            .expect("mixed addition should succeed");
        let normalized_mixed_sum = curve.mixed_add_projective(
            &left_projective.normalize().expect("valid projective point should normalize"),
            &right_affine,
        ).expect("mixed addition should succeed");

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
