use num_bigint::BigUint;

use crate::numerics::cornacchia::{
    CornacchiaError, CornacchiaSolution, cornacchia_candidate_solutions,
    cornacchia_primitive_solutions, cornacchia_with_root,
};

fn bu(value: u64) -> BigUint {
    BigUint::from(value)
}

#[test]
fn cornacchia_finds_a_sum_of_two_squares_from_a_root() {
    let solution =
        cornacchia_with_root(&bu(1), &bu(13), &bu(5)).expect("5² ≡ -1 mod 13 should be valid");
    let solution = solution.expect("13 = 3² + 1·2²");

    assert_eq!(solution.x(), &bu(3));
    assert_eq!(solution.y(), &bu(2));
    assert!(solution.is_primitive());
}

#[test]
fn cornacchia_normalizes_roots_to_the_lower_half_interval() {
    let lower_root =
        cornacchia_with_root(&bu(1), &bu(13), &bu(5)).expect("5 is a valid root modulo 13");
    let upper_root =
        cornacchia_with_root(&bu(1), &bu(13), &bu(8)).expect("8 is a valid root modulo 13");

    assert_eq!(upper_root, lower_root);
}

#[test]
fn cornacchia_returns_none_when_the_final_candidate_is_not_a_square() {
    let solution =
        cornacchia_with_root(&bu(3), &bu(6), &bu(3)).expect("3² ≡ -3 mod 6 should be valid");

    assert_eq!(solution, None);
}

#[test]
fn cornacchia_reports_invalid_inputs_and_roots() {
    assert_eq!(
        cornacchia_with_root(&bu(0), &bu(13), &bu(5)),
        Err(CornacchiaError::ZeroCoefficient)
    );
    assert_eq!(
        cornacchia_with_root(&bu(1), &bu(0), &bu(5)),
        Err(CornacchiaError::ZeroTarget)
    );
    assert_eq!(
        cornacchia_with_root(&bu(1), &bu(1), &bu(0)),
        Err(CornacchiaError::TrivialTarget)
    );
    assert_eq!(
        cornacchia_with_root(&bu(1), &bu(13), &bu(2)),
        Err(CornacchiaError::RootDoesNotSolveCongruence)
    );
}

#[test]
fn candidate_solutions_run_cornacchia_over_all_modular_roots() {
    let solutions =
        cornacchia_candidate_solutions(&bu(1), &bu(65)).expect("-1 has roots modulo 65");

    assert_eq!(
        solutions,
        vec![
            CornacchiaSolution::new(bu(7), bu(4)),
            CornacchiaSolution::new(bu(8), bu(1)),
        ]
    );
}

#[test]
fn candidate_solutions_return_empty_when_no_modular_root_exists() {
    assert_eq!(
        cornacchia_candidate_solutions(&bu(1), &bu(7)),
        Ok(Vec::new())
    );
}

#[test]
fn primitive_solutions_filter_non_primitive_candidates() {
    let candidates = cornacchia_candidate_solutions(&bu(4), &bu(20))
        .expect("-4 has roots modulo 20 and yields two candidates");
    let primitive = cornacchia_primitive_solutions(&bu(4), &bu(20))
        .expect("primitive route should reuse the candidate route");

    assert_eq!(
        candidates,
        vec![
            CornacchiaSolution::new(bu(2), bu(2)),
            CornacchiaSolution::new(bu(4), bu(1)),
        ]
    );
    assert_eq!(primitive, vec![CornacchiaSolution::new(bu(4), bu(1))]);
}
