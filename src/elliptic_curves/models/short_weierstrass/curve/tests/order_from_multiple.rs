use num_bigint::BigUint;

use crate::elliptic_curves::CurveError;

use super::shared::{bu, f7_curve, f7_point};

#[test]
fn point_order_from_multiple_recovers_the_exact_order_by_prime_peeling() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    let report = curve
        .point_order_from_multiple(&point, bu(6), &[(bu(2), 1), (bu(3), 1)])
        .expect("valid annihilating multiple should recover the exact order");

    assert_eq!(report.supplied_multiple(), &bu(6));
    assert_eq!(report.exact_order(), &bu(6));
    assert_eq!(report.remaining_multiple(), &bu(6));
    assert_eq!(report.steps().len(), 2);
    assert_eq!(report.steps()[0].prime(), &bu(2));
    assert_eq!(report.steps()[0].removed_exponent(), 0);
    assert_eq!(report.steps()[1].prime(), &bu(3));
    assert_eq!(report.steps()[1].removed_exponent(), 0);
}

#[test]
fn point_order_from_multiple_can_remove_extra_prime_powers() {
    let curve = f7_curve();
    let point = f7_point(6, 0);

    let report = curve
        .point_order_from_multiple(&point, bu(6), &[(bu(2), 1), (bu(3), 1)])
        .expect("valid annihilating multiple should reduce to the exact order");

    assert_eq!(report.exact_order(), &bu(2));
    assert_eq!(report.steps()[0].removed_exponent(), 0);
    assert_eq!(report.steps()[1].removed_exponent(), 1);
    assert_eq!(report.steps()[1].remaining_multiple_after_step(), &bu(2));
}

#[test]
fn point_order_from_multiple_rejects_invalid_factorizations() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert_eq!(
        curve.point_order_from_multiple(&point, bu(6), &[(bu(4), 1), (bu(3), 1)]),
        Err(CurveError::InvalidPointOrderMultipleFactorization { multiple: bu(6) })
    );
    assert_eq!(
        curve.point_order_from_multiple(&point, bu(6), &[(bu(2), 2)]),
        Err(CurveError::InvalidPointOrderMultipleFactorization { multiple: bu(6) })
    );
}

#[test]
fn point_order_from_multiple_rejects_non_annihilating_multiple() {
    let curve = f7_curve();
    let point = f7_point(2, 1);

    assert_eq!(
        curve.point_order_from_multiple(&point, bu(3), &[(bu(3), 1)]),
        Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint { multiple: bu(3) })
    );
}

#[test]
fn optimized_order_from_multiple_matches_the_baseline_prime_peeling_report() {
    let curve = f7_curve();
    let point = f7_point(6, 0);
    let multiple =
        BigUint::from(2u8).pow(10) * BigUint::from(3u8).pow(8) * BigUint::from(5u8).pow(6);
    let factorization = vec![(bu(2), 10), (bu(3), 8), (bu(5), 6)];

    let optimized = curve
        .point_order_from_multiple(&point, multiple.clone(), &factorization)
        .expect("optimized report should build");
    let baseline = crate::elliptic_curves::models::short_weierstrass::point_order::point_order_from_multiple_baseline(
        &curve,
        &point,
        multiple.clone(),
        &factorization,
    )
    .expect("baseline report should build");

    assert_eq!(optimized, baseline);
}

#[test]
fn trusted_factorization_route_matches_the_validated_route_on_certified_input() {
    let curve = f7_curve();
    let point = f7_point(6, 0);
    let multiple =
        BigUint::from(2u8).pow(10) * BigUint::from(3u8).pow(8) * BigUint::from(5u8).pow(6);
    let factorization = vec![(bu(2), 10), (bu(3), 8), (bu(5), 6)];

    let validated = curve
        .point_order_from_multiple(&point, multiple.clone(), &factorization)
        .expect("validated report should build");
    let trusted = curve
        .point_order_from_multiple_with_trusted_factorization(&point, multiple, &factorization)
        .expect("trusted report should build");

    assert_eq!(trusted, validated);
}
