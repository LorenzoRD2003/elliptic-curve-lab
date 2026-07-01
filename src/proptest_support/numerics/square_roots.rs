use num_bigint::{BigInt, BigUint};
use proptest::prelude::*;

use crate::proptest_support::config::NumericsStrategyConfig;

/// One brute-force-checkable square-root problem modulo `m ≥ 2`.
#[derive(Clone, Debug)]
pub(crate) struct ModularSquareRootCase {
    value: BigInt,
    modulus: BigUint,
    expected_roots: Vec<BigUint>,
}

impl ModularSquareRootCase {
    fn from_i64(value: i64, modulus: u64) -> Self {
        Self {
            value: BigInt::from(value),
            modulus: BigUint::from(modulus),
            expected_roots: brute_force_square_roots_mod_u64(value, modulus),
        }
    }

    pub(crate) fn value(&self) -> &BigInt {
        &self.value
    }

    pub(crate) fn modulus(&self) -> &BigUint {
        &self.modulus
    }

    pub(crate) fn expected_roots(&self) -> &[BigUint] {
        &self.expected_roots
    }
}

/// Generates varied, brute-force-checkable modular square-root cases.
///
/// The sampled moduli deliberately mix primes, prime powers, powers of `2`,
/// squareful composites, and products with several coprime factors. This keeps
/// the Hensel/CRT route honest without making the oracle expensive.
pub(crate) fn arb_modular_square_root_case(
    config: NumericsStrategyConfig,
) -> BoxedStrategy<ModularSquareRootCase> {
    arb_square_root_modulus(config)
        .prop_flat_map(move |modulus| {
            let bound = i64::try_from(config.value_window_factor * modulus)
                .expect("proptest numerics bounds should fit in i64");
            (-bound..=bound).prop_map(move |value| ModularSquareRootCase::from_i64(value, modulus))
        })
        .boxed()
}

pub(crate) fn brute_force_square_roots_mod_u64(value: i64, modulus: u64) -> Vec<BigUint> {
    let target = value.rem_euclid(modulus as i64) as u64;

    (0..modulus)
        .filter(|candidate| (candidate * candidate) % modulus == target)
        .map(BigUint::from)
        .collect()
}

fn arb_square_root_modulus(config: NumericsStrategyConfig) -> BoxedStrategy<u64> {
    let mut moduli = Vec::new();
    moduli.extend(PRIMES);
    moduli.extend(PRIME_POWERS);
    moduli.extend(TWO_POWERS);
    moduli.extend(SQUAREFUL_COMPOSITES);
    moduli.extend(MIXED_COPRIME_PRODUCTS);
    moduli.extend(MANY_FACTOR_COMPOSITES);
    moduli.sort_unstable();
    moduli.dedup();
    moduli.retain(|modulus| *modulus >= 2 && *modulus <= config.max_bruteforce_modulus);

    prop::sample::select(moduli).boxed()
}

const PRIMES: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];

const TWO_POWERS: &[u64] = &[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];

const PRIME_POWERS: &[u64] = &[
    3, 5, 7, 9, 11, 13, 17, 19, 23, 25, 27, 29, 31, 49, 81, 121, 125, 169, 243, 343, 625, 729,
];

const SQUAREFUL_COMPOSITES: &[u64] = &[
    12, 18, 20, 24, 28, 36, 40, 44, 45, 48, 50, 52, 63, 72, 75, 80, 98, 108, 112, 144, 162, 200,
    225, 242, 245, 288, 324, 392, 405, 500, 675, 784, 800, 972, 1125, 1458, 1800, 2025,
];

const MIXED_COPRIME_PRODUCTS: &[u64] = &[
    6, 10, 14, 15, 21, 22, 26, 30, 33, 34, 35, 38, 39, 42, 46, 55, 57, 58, 62, 65, 66, 70, 78, 82,
    85, 91, 95, 102, 105, 110, 114, 130, 154, 165, 182, 195, 210, 231, 330, 390, 462,
];

const MANY_FACTOR_COMPOSITES: &[u64] = &[
    60, 84, 90, 120, 126, 132, 150, 156, 168, 180, 198, 210, 220, 234, 252, 270, 280, 300, 330,
    360, 420, 462, 504, 540, 600, 630, 660, 780, 840, 900, 990, 1020, 1260, 1320, 1560, 1680, 1980,
    2100, 2310,
];
