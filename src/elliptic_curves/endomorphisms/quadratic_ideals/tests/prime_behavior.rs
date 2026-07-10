use crate::elliptic_curves::endomorphisms::{
    quadratic_ideals::{QuadraticPrimeBehavior, QuadraticPrimeBehaviorError},
    quadratic_orders::{ImaginaryQuadraticOrder, QuadraticDiscriminant},
};
use crate::numerics::{PositivePrimeError, positive_mod_biguint};

use super::{bu, maximal_order};

#[test]
fn split_prime_returns_the_two_square_roots_of_the_discriminant() {
    let order = maximal_order(-23);

    let behavior = order
        .prime_behavior(&bu(3))
        .expect("3 is a supported odd prime");

    assert_eq!(
        behavior,
        QuadraticPrimeBehavior::Split {
            roots: (bu(1), bu(2))
        }
    );
}

#[test]
fn split_roots_are_actual_roots_modulo_ell() {
    let order = maximal_order(-23);
    let ell = bu(13);
    let behavior = order
        .prime_behavior(&ell)
        .expect("13 is a supported odd prime");

    let QuadraticPrimeBehavior::Split { roots } = behavior else {
        panic!("13 should split in the order of discriminant -23");
    };
    let discriminant_mod_ell = positive_mod_biguint(order.discriminant().value(), &ell);

    assert!(roots.0 < roots.1);
    assert_eq!(&roots.0 + &roots.1, ell);
    for root in [roots.0, roots.1] {
        assert_eq!((&root * &root) % &ell, discriminant_mod_ell);
    }
}

#[test]
fn inert_prime_reports_no_local_square_root() {
    let order = maximal_order(-23);

    let behavior = order
        .prime_behavior(&bu(5))
        .expect("5 is a supported odd prime");

    assert_eq!(behavior, QuadraticPrimeBehavior::Inert);
}

#[test]
fn ramified_prime_reports_the_repeated_zero_root() {
    let order = maximal_order(-23);

    let behavior = order
        .prime_behavior(&bu(23))
        .expect("23 is a supported odd prime");

    assert_eq!(behavior, QuadraticPrimeBehavior::Ramified { root: bu(0) });
}

#[test]
fn prime_dividing_the_conductor_is_not_treated_as_invertible() {
    let order = ImaginaryQuadraticOrder::new(QuadraticDiscriminant::new(-23), bu(3))
        .expect("conductor 3 should define a non-maximal imaginary quadratic order");

    let behavior = order
        .prime_behavior(&bu(3))
        .expect("3 is prime even though it divides the conductor");

    assert_eq!(
        behavior,
        QuadraticPrimeBehavior::NonInvertibleBecauseDividesConductor
    );
}

#[test]
fn unsupported_dyadic_case_is_explicit_for_now() {
    let order = maximal_order(-23);

    let error = order
        .prime_behavior(&bu(2))
        .expect_err("2 ∤ f should be outside the first odd-prime API");

    assert_eq!(error, QuadraticPrimeBehaviorError::UnsupportedPrimeTwo);
}

#[test]
fn composite_input_is_rejected_before_local_classification() {
    let order = maximal_order(-23);

    let error = order.prime_behavior(&bu(9)).expect_err("9 is not a prime");

    assert_eq!(
        error,
        QuadraticPrimeBehaviorError::InvalidPrime(PositivePrimeError::Composite)
    );
}
