use super::search_bsgs::{
    HasseBsgsConfig, HasseBsgsParity, HasseBsgsTraversal,
    find_annihilating_multiple_in_interval_bsgs,
    find_annihilating_multiple_in_interval_bsgs_with_config,
};
use crate::elliptic_curves::traits::HasseMultipleSearchCurveModel;
use crate::elliptic_curves::{
    CurveModel, EnumerableCurveModel, GroupCurveModel, ShortWeierstrassCurve,
};
use crate::fields::{Field, Fp};

type F241 = Fp<241>;

#[test]
fn bsgs_hasse_search_finds_an_annihilating_multiple_inside_the_same_hasse_interval() {
    let curve =
        ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
            .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval = crate::elliptic_curves::HasseInterval::for_q(241).expect("valid Hasse interval");

    let naive = curve
        .find_annihilating_multiple_in_interval_naive(&point, interval.clone())
        .expect("naive Hasse search should succeed");
    let bsgs = find_annihilating_multiple_in_interval_bsgs(&curve, &point, interval)
        .expect("BSGS Hasse search should succeed")
        .expect("Hasse's theorem guarantees an annihilating multiple");

    assert!(naive.interval().contains(bsgs));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(bsgs).expect("small-prime Hasse candidates fit in u64")
    ));
}

#[test]
fn configurable_bsgs_defaults_preserve_the_current_search_result() {
    let curve =
        ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
            .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval = crate::elliptic_curves::HasseInterval::for_q(241).expect("valid Hasse interval");

    let default_result =
        find_annihilating_multiple_in_interval_bsgs(&curve, &point, interval.clone())
            .expect("default BSGS should succeed");
    let explicit_default = find_annihilating_multiple_in_interval_bsgs_with_config(
        &curve,
        &point,
        interval,
        HasseBsgsConfig {
            traversal: HasseBsgsTraversal::LeftToRight,
            use_fast_negation: true,
            known_parity: HasseBsgsParity::Unknown,
        },
    )
    .expect("configurable BSGS should succeed");

    assert_eq!(default_result, explicit_default);
}

#[test]
fn fast_negation_and_plain_bsgs_find_valid_annihilating_multiples() {
    let curve =
        ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
            .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval = crate::elliptic_curves::HasseInterval::for_q(241).expect("valid Hasse interval");

    let plain = find_annihilating_multiple_in_interval_bsgs_with_config(
        &curve,
        &point,
        interval.clone(),
        HasseBsgsConfig {
            traversal: HasseBsgsTraversal::LeftToRight,
            use_fast_negation: false,
            known_parity: HasseBsgsParity::Unknown,
        },
    )
    .expect("plain BSGS should succeed")
    .expect("plain BSGS should find an annihilating multiple");
    let fast = find_annihilating_multiple_in_interval_bsgs_with_config(
        &curve,
        &point,
        interval.clone(),
        HasseBsgsConfig {
            traversal: HasseBsgsTraversal::LeftToRight,
            use_fast_negation: true,
            known_parity: HasseBsgsParity::Unknown,
        },
    )
    .expect("fast-negation BSGS should succeed")
    .expect("fast-negation BSGS should find an annihilating multiple");

    assert!(interval.contains(plain));
    assert!(interval.contains(fast));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(plain).expect("small-prime Hasse candidates fit in u64")
    ));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(fast).expect("small-prime Hasse candidates fit in u64")
    ));
}
