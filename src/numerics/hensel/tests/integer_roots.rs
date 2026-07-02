use crate::numerics::hensel::{
    HenselIntegerRootSearchConfig, HenselIntegerRootSearchReport, HenselIntegerRootTrace,
    HenselLiftError, find_integer_roots_by_hensel, hensel_lift_integer_root,
};

use super::{bi, bu, integer_polynomial};

#[test]
fn integer_root_lift_certifies_a_positive_root() {
    let polynomial = integer_polynomial(&[-3, -2, 1]);
    let trace: HenselIntegerRootTrace =
        hensel_lift_integer_root(&polynomial, &bi(3), &bu(5), &bu(4))
            .expect("x = 3 is a simple root modulo 5 and an integer root");

    assert_eq!(trace.candidate_root(), &bi(3));
    assert_eq!(trace.root_bound(), &bu(4));
    assert_eq!(trace.modulus(), &bu(25));
    assert_eq!(trace.lift_trace().reached_level(), 2);
    assert_eq!(polynomial.evaluate(trace.candidate_root()), bi(0));
}

#[test]
fn integer_root_lift_certifies_a_negative_root() {
    let polynomial = integer_polynomial(&[-2, -1, 1]);
    let trace = hensel_lift_integer_root(&polynomial, &bi(4), &bu(5), &bu(3))
        .expect("x = -1 is represented by 4 modulo 5 and should be recovered");

    assert_eq!(trace.candidate_root(), &bi(-1));
    assert_eq!(trace.modulus(), &bu(25));
    assert_eq!(polynomial.evaluate(trace.candidate_root()), bi(0));
}

#[test]
fn target_level_one_can_certify_when_prime_already_exceeds_the_bound() {
    let polynomial = integer_polynomial(&[-6, -1, 1]);
    let trace = hensel_lift_integer_root(&polynomial, &bi(3), &bu(11), &bu(3))
        .expect("the initial residue already determines the bounded root");

    assert_eq!(trace.candidate_root(), &bi(3));
    assert_eq!(trace.modulus(), &bu(11));
    assert!(trace.lift_trace().steps().is_empty());
}

#[test]
fn target_precision_is_chosen_from_the_root_bound() {
    let polynomial = integer_polynomial(&[-5, 1]);
    let trace = hensel_lift_integer_root(&polynomial, &bi(1), &bu(2), &bu(10))
        .expect("p^e should grow until it exceeds twice the bound");

    assert_eq!(trace.candidate_root(), &bi(5));
    assert_eq!(trace.modulus(), &bu(32));
    assert_eq!(trace.lift_trace().reached_level(), 5);
}

#[test]
fn centered_representative_recovers_the_negative_integer_lift() {
    let polynomial = integer_polynomial(&[-2, -1, 1]);
    let trace = hensel_lift_integer_root(&polynomial, &bi(4), &bu(5), &bu(3))
        .expect("the residue 4 modulo 5 should lift to the bounded root -1");

    assert_eq!(trace.candidate_root(), &bi(-1));
}

#[test]
fn integer_root_bound_is_inclusive() {
    let polynomial = integer_polynomial(&[-7, 1]);
    let trace = hensel_lift_integer_root(&polynomial, &bi(7), &bu(11), &bu(7))
        .expect("a root exactly on the absolute bound should be accepted");

    assert_eq!(trace.candidate_root(), &bi(7));
    assert_eq!(trace.root_bound(), &bu(7));
}

#[test]
fn non_root_seed_is_rejected_before_integer_recovery() {
    let polynomial = integer_polynomial(&[-3, -2, 1]);

    assert_eq!(
        hensel_lift_integer_root(&polynomial, &bi(2), &bu(5), &bu(4)),
        Err(HenselLiftError::RootDoesNotSolveCurrentModulus)
    );
}

#[test]
fn singular_derivative_seed_is_rejected_for_the_simple_route() {
    let polynomial = integer_polynomial(&[0, 0, 1]);

    assert_eq!(
        hensel_lift_integer_root(&polynomial, &bi(0), &bu(5), &bu(10)),
        Err(HenselLiftError::SingularDerivativeModPrime)
    );
}

#[test]
fn too_small_bound_does_not_certify_the_lifted_integer() {
    let polynomial = integer_polynomial(&[-3, -2, 1]);

    assert_eq!(
        hensel_lift_integer_root(&polynomial, &bi(3), &bu(7), &bu(2)),
        Err(HenselLiftError::IntegerRootNotCertifiedInBound)
    );
}

#[test]
fn matching_residue_that_is_not_an_integer_root_is_not_certified() {
    let polynomial = integer_polynomial(&[-6, 0, 1]);

    assert_eq!(
        hensel_lift_integer_root(&polynomial, &bi(1), &bu(5), &bu(2)),
        Err(HenselLiftError::IntegerRootNotCertifiedInBound)
    );
}

#[test]
fn search_certifies_all_simple_roots_inside_the_bound() {
    let polynomial = integer_polynomial(&[-3, -2, 1]);
    let config = HenselIntegerRootSearchConfig::new(bu(5), bu(4));
    let report: HenselIntegerRootSearchReport = find_integer_roots_by_hensel(&polynomial, config)
        .expect("both bounded integer roots are simple modulo 5");

    assert!(report.has_certified_roots());
    assert_eq!(report.config().prime(), &bu(5));
    assert_eq!(report.config().root_bound(), &bu(4));
    assert_eq!(report.config().max_seed_scan(), 10_000);
    assert_eq!(report.simple_seed_count(), 2);
    assert_eq!(report.singular_seed_count(), 0);
    assert_eq!(report.uncertified_seed_count(), 0);
    assert_eq!(report.certified_roots(), &[bi(-1), bi(3)]);
    assert_eq!(report.traces().len(), 2);
}

#[test]
fn search_can_report_no_simple_roots_for_the_chosen_prime() {
    let polynomial = integer_polynomial(&[1, 0, 1]);
    let config = HenselIntegerRootSearchConfig::new(bu(3), bu(10));
    let report =
        find_integer_roots_by_hensel(&polynomial, config).expect("x^2 + 1 has no roots modulo 3");

    assert!(!report.has_certified_roots());
    assert_eq!(report.simple_seed_count(), 0);
    assert_eq!(report.singular_seed_count(), 0);
    assert_eq!(report.uncertified_seed_count(), 0);
    assert!(report.certified_roots().is_empty());
    assert!(report.traces().is_empty());
}

#[test]
fn search_counts_singular_modular_roots_without_lifting_them() {
    let polynomial = integer_polynomial(&[0, 0, 1]);
    let config = HenselIntegerRootSearchConfig::new(bu(5), bu(10));
    let report = find_integer_roots_by_hensel(&polynomial, config)
        .expect("singular roots are skipped by the simple-root search");

    assert_eq!(report.simple_seed_count(), 0);
    assert_eq!(report.singular_seed_count(), 1);
    assert!(report.certified_roots().is_empty());
}

#[test]
fn search_counts_simple_seeds_that_do_not_certify_integer_roots() {
    let polynomial = integer_polynomial(&[-6, 0, 1]);
    let config = HenselIntegerRootSearchConfig::new(bu(5), bu(2));
    let report = find_integer_roots_by_hensel(&polynomial, config)
        .expect("simple modular roots can fail integer certification");

    assert_eq!(report.simple_seed_count(), 2);
    assert_eq!(report.singular_seed_count(), 0);
    assert_eq!(report.uncertified_seed_count(), 2);
    assert!(report.certified_roots().is_empty());
}

#[test]
fn search_respects_the_explicit_seed_scan_limit() {
    let polynomial = integer_polynomial(&[-3, -2, 1]);
    let config = HenselIntegerRootSearchConfig::new(bu(11), bu(4)).with_max_seed_scan(5);

    assert_eq!(
        find_integer_roots_by_hensel(&polynomial, config),
        Err(HenselLiftError::SeedScanLimitExceeded)
    );
}
