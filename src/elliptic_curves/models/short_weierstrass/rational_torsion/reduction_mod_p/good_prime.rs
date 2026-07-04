use num_prime::nt_funcs::next_prime;

use crate::elliptic_curves::short_weierstrass::rational_torsion::{
    integral_model::{RationalIntegralModel, integral_rational_to_bigint},
    reduction_mod_p::small_prime_field::{ReductionPrime, ReductionResidue},
};

const FIRST_GOOD_REDUCTION_PRIME_CANDIDATE: u32 = 11;

/// A runtime prime where the integral model has good reduction.
///
/// For an integral short-Weierstrass model `E: y² = x³ + Ax + B`, good
/// reduction at `p` means `Δ(E) ≠ 0 mod p`. The stage-5 route records both the
/// chosen prime and this non-zero discriminant residue as the local certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct GoodReductionPrime {
    prime: ReductionPrime,
    discriminant_mod_p: ReductionResidue,
}

impl GoodReductionPrime {
    /// Finds the first prime `p ≥ 11` that does not divide `Δ(E)`.
    ///
    /// Complexity: linear in the number of primes tested, with one `Δ mod p`
    /// reduction per tested prime.
    pub(super) fn first_for_integral_model(model: &RationalIntegralModel) -> Option<Self> {
        let discriminant = integral_rational_to_bigint(&model.curve().discriminant())
            .expect("RationalIntegralModel should have integral discriminant");
        let mut previous = FIRST_GOOD_REDUCTION_PRIME_CANDIDATE - 1;

        // `next_prime(n)` returns the first prime strictly greater than `n`.
        while let Some(candidate) = next_prime(&previous, None) {
            let prime = ReductionPrime::new(candidate).expect("next_prime should return a prime");
            let discriminant_mod_p = prime.reduce_bigint(&discriminant);
            if !discriminant_mod_p.is_zero() {
                return Some(Self {
                    prime,
                    discriminant_mod_p,
                });
            }
            previous = candidate;
        }
        None
    }

    pub(super) fn prime(&self) -> ReductionPrime {
        self.prime
    }

    pub(super) fn discriminant_mod_p(&self) -> ReductionResidue {
        self.discriminant_mod_p
    }
}
