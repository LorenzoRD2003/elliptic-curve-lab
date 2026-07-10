use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;

use super::{bu, maximal_order};

#[test]
fn prime_norm_ideal_records_order_norm_and_split_root() {
    let order = maximal_order(-23);

    let ideal = PrimeNormIdeal::split(order.clone(), bu(3), bu(1))
        .expect("root 1 should select a split prime ideal above 3");

    assert_eq!(ideal.order(), &order);
    assert_eq!(ideal.norm(), &bu(3));
    assert_eq!(ideal.root_mod_ell(), &bu(1));
    assert!(ideal.is_split());
    assert!(!ideal.is_ramified());
}

#[test]
fn prime_norm_ideal_conjugation_switches_to_the_other_root() {
    let order = maximal_order(-23);
    let ideal = PrimeNormIdeal::split(order.clone(), bu(13), bu(4))
        .expect("root 4 should select a split prime ideal above 13");

    let conjugate = ideal.conjugate();

    assert_eq!(conjugate.order(), &order);
    assert_eq!(conjugate.norm(), &bu(13));
    assert_eq!(conjugate.root_mod_ell(), &bu(9));
    assert_eq!(conjugate.conjugate(), ideal);
}

#[test]
fn prime_norm_ideal_wrapper_delegates_basic_data_and_conjugation() {
    let order = maximal_order(-23);
    let ideal = PrimeNormIdeal::split(order.clone(), bu(3), bu(1))
        .expect("root 1 should select a split prime ideal above 3");

    assert_eq!(ideal.order(), &order);
    assert_eq!(ideal.norm(), &bu(3));
    assert_eq!(ideal.conjugate().order(), &order);
    assert_eq!(ideal.conjugate().norm(), &bu(3));
    assert_eq!(ideal.conjugate().root_mod_ell(), &bu(2));
}

#[test]
fn prime_norm_ideal_builds_ramified_prime_with_repeated_root() {
    let order = maximal_order(-23);

    let ideal = PrimeNormIdeal::ramified(order.clone(), bu(23))
        .expect("23 ramifies in the order of discriminant -23");

    assert_eq!(ideal.order(), &order);
    assert_eq!(ideal.norm(), &bu(23));
    assert_eq!(ideal.root_mod_ell(), &bu(0));
    assert!(!ideal.is_split());
    assert!(ideal.is_ramified());
}

#[test]
fn ramified_prime_ideal_is_fixed_by_conjugation() {
    let order = maximal_order(-23);
    let ideal = PrimeNormIdeal::ramified(order, bu(23))
        .expect("23 ramifies in the order of discriminant -23");

    assert_eq!(ideal.conjugate(), ideal);
}
