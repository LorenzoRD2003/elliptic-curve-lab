use crate::elliptic_curves::{
    CurveError,
    frobenius::group_order::{GroupOrderReport, GroupOrderStrategy},
};

use super::shared::F43;
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::{Fp, traits::Field};

#[test]
fn schoof_strategy_reports_when_the_current_route_blocks() {
    let curve = ShortWeierstrassCurve::<crate::fields::Fp<7>>::new(
        crate::fields::Fp::<7>::from_i64(2),
        crate::fields::Fp::<7>::from_i64(3),
    )
    .expect("valid curve");

    assert_eq!(
        curve.group_order_by(GroupOrderStrategy::Schoof),
        Err(CurveError::SchoofBlockedOnOddPrime { odd_prime: 3 })
    );
}

#[test]
fn schoof_strategy_matches_exhaustive_on_one_f241_curve() {
    type F241 = Fp<241>;

    for a in -8..=8 {
        for b in -8..=8 {
            let Ok(curve) =
                ShortWeierstrassCurve::<F241>::new(F241::from_i64(a), F241::from_i64(b))
            else {
                continue;
            };
            let Ok(schoof) = curve.group_order_by(GroupOrderStrategy::Schoof) else {
                continue;
            };
            let exhaustive = curve
                .group_order_by(GroupOrderStrategy::Exhaustive)
                .expect("exhaustive route should succeed over F241");

            assert_eq!(schoof.curve_order(), exhaustive.curve_order());
            assert_eq!(schoof.trace(), exhaustive.trace());

            let GroupOrderReport::Schoof(report) = schoof else {
                panic!("Schoof strategy should preserve its own report variant");
            };

            assert!(!report.attempted_odd_primes().is_empty());
            assert!(*report.combined_crt_modulus() > num_bigint::BigUint::from(0u8));
            return;
        }
    }

    panic!("expected to find one small F241 curve resolved by the current Schoof strategy");
}

#[test]
fn schoof_strategy_can_be_found_on_a_small_f43_curve() {
    for a in -3..=3 {
        for b in -3..=3 {
            let Ok(curve) = ShortWeierstrassCurve::<F43>::new(F43::from_i64(a), F43::from_i64(b))
            else {
                continue;
            };
            let Ok(schoof) = curve.group_order_by(GroupOrderStrategy::Schoof) else {
                continue;
            };
            let exhaustive = curve
                .group_order_by(GroupOrderStrategy::Exhaustive)
                .expect("exhaustive route should succeed over F43");

            assert_eq!(schoof.curve_order(), exhaustive.curve_order());
            assert_eq!(schoof.trace(), exhaustive.trace());
            return;
        }
    }

    panic!("expected at least one small F43 curve to resolve under the current Schoof strategy");
}
