use num_bigint::BigInt;
use num_rational::BigRational;
use proptest::prelude::*;

use crate::proptest_support::config::FieldStrategyConfig;

/// Returns a strategy for exact rational numbers with bounded numerator and
/// positive denominator.
pub fn arb_q_elem(config: FieldStrategyConfig) -> BoxedStrategy<BigRational> {
    let max_abs = config.max_abs_i64.max(1);
    (-max_abs..=max_abs, 1..=max_abs)
        .prop_map(|(numerator, denominator)| {
            BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
        })
        .boxed()
}
