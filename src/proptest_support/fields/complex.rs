use num_complex::Complex64;
use proptest::prelude::*;

use crate::proptest_support::config::FieldStrategyConfig;

/// Returns a bounded approximate complex sample.
pub fn arb_complex_approx(config: FieldStrategyConfig) -> BoxedStrategy<Complex64> {
    let max_real = config.max_real_norm;
    let max_imag = config.max_imaginary_norm;

    (-max_real..=max_real, -max_imag..=max_imag)
        .prop_map(|(re, im)| Complex64::new(re, im))
        .boxed()
}
