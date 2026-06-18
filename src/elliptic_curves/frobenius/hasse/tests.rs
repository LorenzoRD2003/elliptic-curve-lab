use crate::elliptic_curves::traits::HasseIntervalSearchCurveModel;
use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::hasse::search::{HasseBsgsConfig, HasseBsgsParity, HasseBsgsTraversal},
    models::short_weierstrass::group_order_parity::GroupOrderParity,
    traits::{AffineCurveModel, CurveModel, EnumerableCurveModel, GroupCurveModel},
};
use crate::fields::{
    Fp,
    traits::{Field, FiniteField},
};

type F241 = Fp<241>;

#[test]
fn hasse_interval_can_be_built_directly_from_a_field_family() {
    let from_field = crate::elliptic_curves::frobenius::HasseInterval::for_field::<F241>()
        .expect("valid Hasse interval");
    let from_order = crate::elliptic_curves::frobenius::HasseInterval::for_q(F241::order())
        .expect("valid Hasse interval");

    assert_eq!(from_field.q(), 241);
    assert_eq!(from_field, from_order);
}

#[test]
fn bsgs_hasse_search_finds_an_annihilating_multiple_inside_the_same_hasse_interval() {
    let curve = ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
        .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval =
        crate::elliptic_curves::frobenius::HasseInterval::for_q(241).expect("valid Hasse interval");

    let naive = curve
        .find_annihilating_multiple_in_interval_naive(&point, interval.clone())
        .expect("naive Hasse search should succeed");
    let bsgs = curve
        .find_annihilating_multiple_in_interval_bsgs(&point, interval)
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
    let curve = ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
        .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval =
        crate::elliptic_curves::frobenius::HasseInterval::for_q(241).expect("valid Hasse interval");

    let default_result =
        HasseIntervalSearchCurveModel::find_annihilating_multiple_in_interval_bsgs(
            &curve,
            &point,
            interval.clone(),
        )
        .expect("default BSGS should succeed");
    let explicit_default = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval,
            HasseBsgsConfig::new()
                .with_traversal(HasseBsgsTraversal::LeftToRight)
                .with_fast_negation(true)
                .with_known_parity(HasseBsgsParity::Unknown),
        )
        .expect("configurable BSGS should succeed");

    assert_eq!(default_result, explicit_default);
}

#[test]
fn fast_negation_and_plain_bsgs_find_valid_annihilating_multiples() {
    let curve = ShortWeierstrassCurve::<F241>::new(F241::from_i64(2), F241::from_i64(3))
        .expect("valid curve");
    let point = curve
        .points()
        .into_iter()
        .find(|point| !curve.is_identity(point))
        .expect("small finite curve should contain a non-identity point");
    let interval =
        crate::elliptic_curves::frobenius::HasseInterval::for_q(241).expect("valid Hasse interval");

    let plain = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::new()
                .with_traversal(HasseBsgsTraversal::LeftToRight)
                .with_fast_negation(false)
                .with_known_parity(HasseBsgsParity::Unknown),
        )
        .expect("plain BSGS should succeed")
        .expect("plain BSGS should find an annihilating multiple");
    let fast = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::new()
                .with_traversal(HasseBsgsTraversal::LeftToRight)
                .with_fast_negation(true)
                .with_known_parity(HasseBsgsParity::Unknown),
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

#[test]
fn known_even_parity_and_unknown_bsgs_find_valid_annihilating_multiples() {
    let curve = ShortWeierstrassCurve::<F241>::new(F241::zero(), F241::one()).expect("valid curve");
    let point = curve
        .point(F241::zero(), F241::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval =
        crate::elliptic_curves::frobenius::HasseInterval::for_q(241).expect("valid Hasse interval");

    assert_eq!(
        curve.group_order_parity_from_two_torsion(),
        GroupOrderParity::Even
    );

    let unknown = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::default(),
        )
        .expect("unknown-parity BSGS should succeed")
        .expect("unknown-parity BSGS should find an annihilating multiple");
    let even = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::new().with_known_parity(HasseBsgsParity::Even),
        )
        .expect("even-parity BSGS should succeed")
        .expect("even-parity BSGS should find an annihilating multiple");

    assert!(interval.contains(unknown));
    assert!(interval.contains(even));
    assert_eq!(even % 2, 0);
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(unknown).expect("small-prime Hasse candidates fit in u64")
    ));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(even).expect("small-prime Hasse candidates fit in u64")
    ));
}

#[test]
fn known_odd_parity_and_unknown_bsgs_find_valid_annihilating_multiples() {
    let curve =
        ShortWeierstrassCurve::<F241>::new(F241::from_i64(3), F241::one()).expect("valid curve");
    let point = curve
        .point(F241::zero(), F241::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval =
        crate::elliptic_curves::frobenius::HasseInterval::for_q(241).expect("valid Hasse interval");

    assert_eq!(
        curve.group_order_parity_from_two_torsion(),
        GroupOrderParity::Odd
    );

    let unknown = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::default(),
        )
        .expect("unknown-parity BSGS should succeed")
        .expect("unknown-parity BSGS should find an annihilating multiple");
    let odd = curve
        .find_annihilating_multiple_in_interval_bsgs_with_config(
            &point,
            interval.clone(),
            HasseBsgsConfig::new().with_known_parity(HasseBsgsParity::Odd),
        )
        .expect("odd-parity BSGS should succeed")
        .expect("odd-parity BSGS should find an annihilating multiple");

    assert!(interval.contains(unknown));
    assert!(interval.contains(odd));
    assert_eq!(odd % 2, 1);
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(unknown).expect("small-prime Hasse candidates fit in u64")
    ));
    assert!(curve.is_torsion_point(
        &point,
        u64::try_from(odd).expect("small-prime Hasse candidates fit in u64")
    ));
}
