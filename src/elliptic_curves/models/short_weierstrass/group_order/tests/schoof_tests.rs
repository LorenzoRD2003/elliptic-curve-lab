use crate::elliptic_curves::frobenius::group_order::{
    FiniteFieldGroupOrderStrategy, GroupOrderReport, SmallFieldGroupOrderStrategy,
};
use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::ShortWeierstrassCurve;

#[test]
fn schoof_strategy_can_skip_one_blocked_prime_and_still_resolve() {
    let curve = ShortWeierstrassCurve::<crate::fields::Fp7>::new(
        crate::fields::Fp7::from_i64(2),
        crate::fields::Fp7::from_i64(3),
    )
    .expect("valid curve");

    let report = curve
        .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
        .expect("the automatic finite-field Schoof route should continue after one blocked prime");

    let GroupOrderReport::Schoof(report) = report else {
        panic!("Schoof strategy should preserve its own report variant");
    };

    assert_eq!(report.resolved().curve_order(), BigUint::from(6u8));
    assert_eq!(report.resolved().trace(), BigInt::from(2));
    assert_eq!(report.attempted_odd_primes(), &[3, 5, 11]);
    assert_eq!(
        report.combined_crt_modulus(),
        &num_bigint::BigUint::from(22u8)
    );
}

#[test]
fn schoof_strategy_matches_exhaustive_on_one_f241_curve() {
    type F241 = crate::fields::Fp241;
    let curve = ShortWeierstrassCurve::<F241>::new(F241::from_i64(-4), F241::from_i64(-4))
        .expect("benchmark F241 curve should be valid");
    let schoof = curve
        .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
        .expect("the benchmark F241 curve should resolve under the current Schoof strategy");
    let exhaustive = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Exhaustive)
        .expect("exhaustive route should succeed over F241");

    assert_eq!(schoof.curve_order(), exhaustive.curve_order());
    assert_eq!(schoof.trace(), exhaustive.trace());

    let GroupOrderReport::Schoof(report) = schoof else {
        panic!("Schoof strategy should preserve its own report variant");
    };

    assert!(!report.attempted_odd_primes().is_empty());
    assert!(*report.combined_crt_modulus() > num_bigint::BigUint::from(0u8));
}

#[test]
fn schoof_strategy_matches_exhaustive_on_one_small_f7_curve_after_skipping_prime_three() {
    let curve = ShortWeierstrassCurve::<crate::fields::Fp7>::new(
        crate::fields::Fp7::from_i64(2),
        crate::fields::Fp7::from_i64(3),
    )
    .expect("valid curve");
    let schoof = curve
        .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
        .expect("the automatic finite-field Schoof route should resolve this curve");
    let exhaustive = curve
        .group_order_by_small_field(SmallFieldGroupOrderStrategy::Exhaustive)
        .expect("exhaustive route should succeed over F7");

    assert_eq!(schoof.curve_order(), exhaustive.curve_order());
    assert_eq!(schoof.trace(), exhaustive.trace());
}
