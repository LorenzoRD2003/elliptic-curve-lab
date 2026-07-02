use num_bigint::{BigInt, BigUint};

/// Bezout data used by the root formula in the nontrivial `r`-Sylow case.
///
/// When `r | |G|`, there exist integers `s, t` such that
/// `s a + t r^(k+1) = 1`. If `α = aγ = xδ` and `r | x`, the returned root
/// is `ρ = s(x/r)δ + tβ`, where `β = r^k γ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CyclicPrimeRootBezout {
    s: BigInt,
    t: BigInt,
    cofactor: BigUint,
    next_sylow_order: BigUint,
}

impl CyclicPrimeRootBezout {
    pub(crate) fn new(s: BigInt, t: BigInt, cofactor: BigUint, next_sylow_order: BigUint) -> Self {
        Self {
            s,
            t,
            cofactor,
            next_sylow_order,
        }
    }

    /// Returns `s` in `s a + t r^(k+1) = 1`.
    pub fn s(&self) -> &BigInt {
        &self.s
    }

    /// Returns `t` in `s a + t r^(k+1) = 1`.
    pub fn t(&self) -> &BigInt {
        &self.t
    }

    /// Returns the factor `a`.
    pub fn cofactor(&self) -> &BigUint {
        &self.cofactor
    }

    /// Returns `r^(k+1)`.
    pub fn next_sylow_order(&self) -> &BigUint {
        &self.next_sylow_order
    }
}
