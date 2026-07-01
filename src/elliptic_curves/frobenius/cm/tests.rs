use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    endomorphisms::quadratic_orders::QuadraticDiscriminant,
    frobenius::cm::{CmTraceCandidateError, CmTraceSignCurveModel, cm_absolute_trace_candidates},
    traits::{AffineCurveModel, CurveModel, EnumerableCurveModel},
};
use crate::fields::{Fp, traits::Field};

type F43 = Fp<43>;

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn discriminant(value: i64) -> QuadraticDiscriminant {
    QuadraticDiscriminant::new(BigInt::from(value))
}

fn f43_curve(a: i64, b: i64) -> ShortWeierstrassCurve<F43> {
    ShortWeierstrassCurve::<F43>::new(F43::from_i64(a), F43::from_i64(b))
        .expect("test curve should be nonsingular")
}

fn candidate_coordinates(
    candidates: &[crate::elliptic_curves::frobenius::cm::CmTraceCandidate],
) -> Vec<(BigUint, BigUint)> {
    candidates
        .iter()
        .map(|candidate| {
            (
                candidate.absolute_trace().clone(),
                candidate.cm_multiplier().clone(),
            )
        })
        .collect()
}

#[test]
fn cm_absolute_trace_candidates_use_cornacchia_candidates_for_four_p() {
    let candidates =
        cm_absolute_trace_candidates(&discriminant(-7), &bu(29)).expect("4·29 = 2² + 7·4²");

    assert_eq!(candidate_coordinates(&candidates), vec![(bu(2), bu(4))]);
}

#[test]
fn cm_absolute_trace_candidates_include_direct_four_p_representations() {
    let candidates =
        cm_absolute_trace_candidates(&discriminant(-11), &bu(3)).expect("4·3 = 1² + 11·1²");

    assert_eq!(candidate_coordinates(&candidates), vec![(bu(1), bu(1))]);
}

#[test]
fn cm_absolute_trace_candidates_return_empty_when_no_representation_exists() {
    let candidates = cm_absolute_trace_candidates(&discriminant(-7), &bu(3))
        .expect("negative discriminant and positive p are valid inputs");

    assert!(candidates.is_empty());
}

#[test]
fn cm_absolute_trace_candidates_reject_invalid_inputs() {
    assert_eq!(
        cm_absolute_trace_candidates(&discriminant(0), &bu(29)),
        Err(CmTraceCandidateError::NonNegativeDiscriminant)
    );
    assert_eq!(
        cm_absolute_trace_candidates(&discriminant(5), &bu(29)),
        Err(CmTraceCandidateError::NonNegativeDiscriminant)
    );
    assert_eq!(
        cm_absolute_trace_candidates(&discriminant(-7), &bu(0)),
        Err(CmTraceCandidateError::ZeroPrime)
    );
}

#[test]
fn cm_trace_sign_test_recovers_positive_trace_from_distinguishing_point() {
    let curve = f43_curve(1, 1);
    let point = curve
        .point(F43::zero(), F43::one())
        .expect("(0, 1) lies on y² = x³ + x + 1 over F43");

    assert_eq!(curve.order(), 34);
    assert_eq!(
        curve.cm_trace_from_absolute_trace_with_point(&bu(43), &bu(10), &point),
        Ok(Some(BigInt::from(10)))
    );
}

#[test]
fn cm_trace_sign_test_recovers_negative_trace_from_distinguishing_point() {
    let curve = f43_curve(2, 3);
    let point = curve
        .point(F43::one(), F43::from_i64(7))
        .expect("(1, 7) lies on y² = x³ + 2x + 3 over F43");

    assert_eq!(curve.order(), 46);
    assert_eq!(
        curve.cm_trace_from_absolute_trace_with_point(&bu(43), &bu(2), &point),
        Ok(Some(BigInt::from(-2)))
    );
}

#[test]
fn cm_trace_sign_test_returns_none_when_point_does_not_distinguish_sign() {
    let curve = f43_curve(1, 1);
    let identity = curve.identity();

    assert_eq!(
        curve.cm_trace_from_absolute_trace_with_point(&bu(43), &bu(10), &identity),
        Ok(None)
    );
}

#[test]
fn cm_trace_sign_test_can_use_enumeration_as_a_deterministic_witness_source() {
    let curve = f43_curve(2, 3);

    assert_eq!(
        curve.cm_trace_from_absolute_trace_by_enumeration(&bu(43), &bu(2)),
        Ok(Some(BigInt::from(-2)))
    );
}

#[test]
fn cm_trace_sign_test_can_use_sampler_backed_random_points() {
    let curve = f43_curve(1, 1);
    let mut indices = [0usize, 1].into_iter();
    let mut sampler = |upper_bound: usize| indices.next().filter(|index| *index < upper_bound);

    assert_eq!(
        curve.cm_trace_from_absolute_trace_by_random_points(&bu(43), &bu(10), &mut sampler, 2),
        Ok(Some(BigInt::from(10)))
    );
}

#[test]
fn cm_trace_sign_test_random_points_honor_attempt_limit() {
    let curve = f43_curve(1, 1);
    let mut sampler = |_| Some(1usize);

    assert_eq!(
        curve.cm_trace_from_absolute_trace_by_random_points(&bu(43), &bu(10), &mut sampler, 0),
        Ok(None)
    );
}
