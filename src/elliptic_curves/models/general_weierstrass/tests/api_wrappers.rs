use crate::fields::traits::*;
use num_bigint::BigUint;

use crate::elliptic_curves::{
    GeneralWeierstrassCurve,
    frobenius::group_order::{
        FiniteFieldGroupOrderStrategy, GroupOrderReport, GroupOrderRoute,
        SmallFieldGroupOrderStrategy,
    },
    short_weierstrass::{
        group_exponent::{GroupExponentReport, GroupExponentStrategy},
        point_order::{PointOrderReport, PointOrderStrategy},
    },
    traits::{AffineCurveModel, CurveModelConversion, EnumerableCurveModel, PointIndexSampler},
};

use super::shared::{F2, F3, F5};

#[derive(Clone, Debug)]
struct FixedIndexSampler {
    indices: Vec<usize>,
    cursor: usize,
}

impl FixedIndexSampler {
    fn new(indices: Vec<usize>) -> Self {
        Self { indices, cursor: 0 }
    }
}

impl PointIndexSampler for FixedIndexSampler {
    fn sample_index(&mut self, upper_bound: usize) -> Option<usize> {
        let index = *self.indices.get(self.cursor)?;
        self.cursor += 1;
        (index < upper_bound).then_some(index)
    }
}

#[test]
fn group_order_by_delegates_to_the_short_companion_when_available() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(
        curve
            .group_order_by(FiniteFieldGroupOrderStrategy::Auto)
            .expect("finite-field wrapper should succeed")
            .curve_order(),
        conversion
            .target()
            .group_order_by(FiniteFieldGroupOrderStrategy::Auto)
            .expect("short companion route should succeed")
            .curve_order(),
    );
}

#[test]
fn frobenius_trace_by_delegates_to_the_short_companion_when_available() {
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(
        curve
            .frobenius_trace_by(FiniteFieldGroupOrderStrategy::Auto)
            .expect("finite-field trace wrapper should succeed"),
        conversion
            .target()
            .frobenius_trace_by(FiniteFieldGroupOrderStrategy::Auto)
            .expect("short companion route should succeed"),
    );
}

#[test]
fn finite_field_wrappers_remain_honest_about_short_reduction_limits() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    assert!(matches!(
        curve.group_order_by(FiniteFieldGroupOrderStrategy::Auto),
        Err(crate::elliptic_curves::CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(2u8)
    ));
    assert!(matches!(
        curve.frobenius_trace_by(FiniteFieldGroupOrderStrategy::Schoof),
        Err(crate::elliptic_curves::CurveError::UnsupportedCharacteristic { characteristic })
            if characteristic == BigUint::from(2u8)
    ));
}

#[test]
fn group_order_by_small_field_exhaustive_is_native_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Exhaustive)
        .expect("native exhaustive small-field route should succeed");

    assert!(matches!(report, GroupOrderReport::ExhaustiveTrace(_)));
    assert_eq!(report.route(), GroupOrderRoute::Exhaustive);
    assert_eq!(report.curve_order(), BigUint::from(2u8));
}

#[test]
fn group_order_by_small_field_auto_falls_back_to_exhaustive_when_short_routes_are_unavailable() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");

    let report = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Auto)
        .expect("auto route should fall back to exhaustive");

    assert!(matches!(report, GroupOrderReport::ExhaustiveTrace(_)));
    assert_eq!(report.route(), GroupOrderRoute::Exhaustive);
}

#[test]
fn group_order_by_small_field_quadratic_character_still_matches_the_short_companion_when_supported()
{
    let curve =
        GeneralWeierstrassCurve::<F5>::new(F5::one(), F5::one(), F5::one(), F5::one(), F5::zero())
            .expect("non-singular curve");
    let conversion = curve
        .conversion_to_short_weierstrass()
        .expect("characteristic five should support the short reduction");

    assert_eq!(
        curve
            .group_order_by_small_field(SmallFieldGroupOrderStrategy::QuadraticCharacter)
            .expect("general wrapper should succeed")
            .curve_order(),
        conversion
            .target()
            .group_order_by_small_field(SmallFieldGroupOrderStrategy::QuadraticCharacter)
            .expect("short companion route should succeed")
            .curve_order(),
    );
}

#[test]
fn point_order_by_exhaustive_is_native_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let point = curve
        .point(F2::one(), F2::zero())
        .expect("sample point should lie on the curve");

    let report = curve
        .point_order_by(&point, PointOrderStrategy::Exhaustive)
        .expect("native exhaustive point-order route should succeed");

    let PointOrderReport::Exhaustive(report) = report else {
        panic!("expected the exhaustive point-order variant");
    };
    assert_eq!(report.exact_order(), &BigUint::from(2u8));
}

#[test]
fn point_order_by_from_known_multiple_is_native_in_characteristic_three() {
    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let point = curve
        .point(F3::zero(), F3::zero())
        .expect("sample point should lie on the curve");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::FromKnownMultiple {
                multiple: BigUint::from(6u8),
                factorization: vec![(BigUint::from(2u8), 1), (BigUint::from(3u8), 1)],
            },
        )
        .expect("native order-from-multiple route should succeed");

    let PointOrderReport::FromKnownMultiple(report) = report else {
        panic!("expected the from-known-multiple variant");
    };
    assert_eq!(report.exact_order(), &BigUint::from(6u8));
}

#[test]
fn point_order_by_hasse_interval_naive_uses_the_requested_small_field_route() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let point = curve
        .point(F2::one(), F2::zero())
        .expect("sample point should lie on the curve");

    let report = curve
        .point_order_by(
            &point,
            PointOrderStrategy::HasseIntervalNaive {
                group_order_strategy: SmallFieldGroupOrderStrategy::Exhaustive,
            },
        )
        .expect("native Hasse-interval route should succeed");

    let PointOrderReport::HasseIntervalNaive(report) = report else {
        panic!("expected the HasseIntervalNaive variant");
    };
    assert_eq!(report.exact_order(), &BigUint::from(2u8));
    assert_eq!(
        report.group_order_report().route(),
        GroupOrderRoute::Exhaustive
    );
}

#[test]
fn group_exponent_by_exhaustive_is_native_in_characteristic_two() {
    let curve =
        GeneralWeierstrassCurve::<F2>::new(F2::one(), F2::zero(), F2::one(), F2::zero(), F2::one())
            .expect("non-singular curve in characteristic two");
    let mut sampler = FixedIndexSampler::new(vec![]);

    let report = curve
        .group_exponent_by(GroupExponentStrategy::Exhaustive, &mut sampler)
        .expect("native exhaustive exponent route should succeed");

    let GroupExponentReport::Exhaustive(exponent) = report else {
        panic!("expected the exhaustive exponent variant");
    };
    assert_eq!(exponent, BigUint::from(2u8));
}

#[test]
fn group_exponent_by_random_points_uses_point_order_by_natively() {
    let curve =
        GeneralWeierstrassCurve::<F3>::new(F3::one(), F3::zero(), F3::one(), F3::one(), F3::zero())
            .expect("non-singular curve in characteristic three");
    let point = curve
        .point(F3::zero(), F3::zero())
        .expect("sample point should lie on the curve");
    let index = curve
        .points()
        .iter()
        .position(|candidate| candidate == &point)
        .expect("sample point should appear in the enumerated point set");
    let mut sampler = FixedIndexSampler::new(vec![index]);

    let report = curve
        .group_exponent_by(
            GroupExponentStrategy::RandomPoints {
                max_samples: 1,
                point_order_strategy: PointOrderStrategy::Exhaustive,
            },
            &mut sampler,
        )
        .expect("native random-point exponent route should succeed");

    let GroupExponentReport::RandomPoints(report) = report else {
        panic!("expected the random-points exponent variant");
    };
    assert_eq!(report.exponent_lower_bound(), &BigUint::from(6u8));
}
