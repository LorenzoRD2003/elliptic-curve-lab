use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::{
    endomorphisms::quadratic_orders::QuadraticDiscriminant,
    frobenius::cm::{CmTraceCandidateError, cm_absolute_trace_candidates},
};

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

fn discriminant(value: i64) -> QuadraticDiscriminant {
    QuadraticDiscriminant::new(BigInt::from(value))
}

#[test]
fn cm_absolute_trace_candidates_use_cornacchia_candidates_for_four_p() {
    let candidates =
        cm_absolute_trace_candidates(&discriminant(-7), &bu(29)).expect("4·29 = 2² + 7·4²");

    assert_eq!(candidates, vec![bu(2)]);
}

#[test]
fn cm_absolute_trace_candidates_include_direct_four_p_representations() {
    let candidates =
        cm_absolute_trace_candidates(&discriminant(-11), &bu(3)).expect("4·3 = 1² + 11·1²");

    assert_eq!(candidates, vec![bu(1)]);
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
