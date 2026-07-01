use crate::numerics::hensel::{
    HenselLiftError, hensel_lift_simple_root, hensel_lift_simple_root_step,
    polynomial::evaluate_polynomial,
};

use super::{bi, bu, polynomial};

#[test]
fn simple_step_lifts_a_square_root_modulo_the_next_prime_power() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let step = hensel_lift_simple_root_step(&coefficients, &bi(3), &bu(7), 1)
        .expect("3^2 = 2 mod 7 should lift uniquely");

    assert_eq!(step.level(), 1);
    assert_eq!(step.root_before(), &bi(3));
    assert_eq!(step.lift_digit(), &bu(1));
    assert_eq!(step.root_after(), &bi(10));
    assert_eq!(
        evaluate_polynomial(&coefficients, step.root_after()) % bi(49),
        bi(0)
    );
}

#[test]
fn repeated_lift_records_each_hensel_digit() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let trace = hensel_lift_simple_root(&coefficients, &bi(3), &bu(7), 4)
        .expect("simple square root should lift through p^4");

    assert_eq!(trace.prime(), &bu(7));
    assert_eq!(trace.coefficients(), coefficients.as_slice());
    assert_eq!(trace.initial_root(), &bi(3));
    assert_eq!(trace.reached_level(), 4);
    assert_eq!(trace.final_root(), &bi(2166));
    assert_eq!(trace.steps().len(), 3);
    assert_eq!(
        trace
            .steps()
            .iter()
            .map(|step| step.lift_digit().clone())
            .collect::<Vec<_>>(),
        vec![bu(1), bu(2), bu(6)]
    );
    assert_eq!(
        evaluate_polynomial(&coefficients, trace.final_root()) % bi(2401),
        bi(0)
    );
}

#[test]
fn target_level_one_validates_and_returns_the_normalized_initial_root() {
    let coefficients = polynomial(&[-2, 0, 1]);
    let trace = hensel_lift_simple_root(&coefficients, &bi(-4), &bu(7), 1)
        .expect("-4 is the same root as 3 modulo 7");

    assert_eq!(trace.initial_root(), &bi(3));
    assert_eq!(trace.final_root(), &bi(3));
    assert_eq!(trace.reached_level(), 1);
    assert!(trace.steps().is_empty());
}

#[test]
fn non_root_is_rejected_before_lifting() {
    let coefficients = polynomial(&[-2, 0, 1]);

    assert_eq!(
        hensel_lift_simple_root(&coefficients, &bi(2), &bu(7), 3),
        Err(HenselLiftError::RootDoesNotSolveCurrentModulus)
    );
}

#[test]
fn singular_derivative_is_rejected_for_the_simple_route() {
    let coefficients = polynomial(&[0, 0, 1]);

    assert_eq!(
        hensel_lift_simple_root(&coefficients, &bi(0), &bu(5), 2),
        Err(HenselLiftError::SingularDerivativeModPrime)
    );
}

#[test]
fn invalid_inputs_are_reported_explicitly() {
    assert_eq!(
        hensel_lift_simple_root(&[], &bi(0), &bu(5), 2),
        Err(HenselLiftError::EmptyPolynomial)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[1]), &bi(0), &bu(5), 2),
        Err(HenselLiftError::ConstantPolynomial)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(1), 2),
        Err(HenselLiftError::NonPrimeModulus)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(9), 2),
        Err(HenselLiftError::NonPrimeModulus)
    );
    assert_eq!(
        hensel_lift_simple_root(&polynomial(&[-2, 0, 1]), &bi(3), &bu(7), 0),
        Err(HenselLiftError::ZeroTargetLevel)
    );
}
