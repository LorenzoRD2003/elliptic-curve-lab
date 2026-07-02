use crate::fields::traits::*;
use crate::proptest_support::{
    config::CurveStrategyConfig, elliptic_curves::arb_nonsingular_curve,
};
use proptest::prelude::*;

use super::shared::F17;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(28))]

    #[test]
    fn property_short_weierstrass_invariants_satisfy_the_classical_relation(
        curve in arb_nonsingular_curve::<crate::fields::Fp17>(CurveStrategyConfig::default()),
    ) {
        let left = F17::sub(&F17::cube(&curve.c4()), &F17::square(&curve.c6()));
        let right = F17::mul(&F17::from_i64(1728), &curve.discriminant());

        prop_assert!(F17::eq(&left, &right));
    }
}
