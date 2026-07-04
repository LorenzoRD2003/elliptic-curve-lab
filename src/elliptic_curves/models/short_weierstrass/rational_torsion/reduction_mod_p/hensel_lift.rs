use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::short_weierstrass::rational_torsion::reduction_mod_p::{
    small_prime_field::{ReductionPrime, ReductionResidue},
    torsion_polynomial::{TorsionXPolynomial, TorsionXPolynomialSource},
};
use crate::numerics::hensel::{HenselIntegerRootTrace, HenselLiftError, hensel_lift_integer_root};
use crate::polynomials::IntegerPolynomial;

/// Outcome of attempting to lift one modular `x`-seed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum TorsionXSeedLiftOutcome {
    /// The seed is a multiple root modulo `p`, so the simple Hensel route does
    /// not apply.
    SingularModuloPrime,
    /// The seed lifted p-adically, but did not certify an integer root inside
    /// the Cauchy bound.
    NotCertifiedInBound,
    /// The seed certified an integer root.
    Certified(HenselIntegerRootTrace),
}

/// Report for one modular seed considered by the Hensel stage.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct TorsionXSeedLiftReport {
    #[cfg(test)]
    seed: ReductionResidue,
    outcome: TorsionXSeedLiftOutcome,
}

/// Hensel-lift report for one torsion `x`-criterion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct TorsionXPolynomialLiftReport {
    #[cfg(test)]
    order: usize,
    #[cfg(test)]
    source: TorsionXPolynomialSource,
    #[cfg(test)]
    prime: ReductionPrime,
    #[cfg(test)]
    root_bound: BigUint,
    seeds: Vec<TorsionXSeedLiftReport>,
}

impl TorsionXPolynomial {
    /// Lifts the modular roots of this `x`-criterion to certified integer roots.
    ///
    /// The root bound is the Cauchy bound of the stored primitive polynomial.
    /// Each root of the criterion in `𝔽_p` is handled independently, so
    /// singular and uncertified seeds are recorded as honest outcomes rather
    /// than fatal errors.
    ///
    /// Complexity: `Θ(p·n)`, where `p` is the reduction prime and `n` is the
    /// dense coefficient count of the `x`-criterion.
    pub(super) fn lift_roots_by_hensel(
        &self,
        prime: ReductionPrime,
    ) -> Result<TorsionXPolynomialLiftReport, HenselLiftError> {
        lift_polynomial_roots_by_hensel(
            self.order(),
            self.source(),
            self.polynomial(),
            prime,
            self.roots_mod_prime(prime),
        )
    }
}

impl TorsionXPolynomialLiftReport {
    #[cfg(test)]
    pub(super) fn order(&self) -> usize {
        self.order
    }

    #[cfg(test)]
    pub(super) fn source(&self) -> TorsionXPolynomialSource {
        self.source
    }

    #[cfg(test)]
    pub(super) fn prime(&self) -> ReductionPrime {
        self.prime
    }

    #[cfg(test)]
    pub(super) fn root_bound(&self) -> &BigUint {
        &self.root_bound
    }

    #[cfg(test)]
    pub(super) fn seeds(&self) -> &[TorsionXSeedLiftReport] {
        &self.seeds
    }

    pub(super) fn certified_roots(&self) -> Vec<BigInt> {
        let mut roots = self
            .seeds
            .iter()
            .filter_map(|seed| match seed.outcome() {
                TorsionXSeedLiftOutcome::Certified(trace) => Some(trace.candidate_root().clone()),
                _ => None,
            })
            .collect::<Vec<_>>();
        roots.sort();
        roots.dedup();
        roots
    }

    #[cfg(test)]
    pub(super) fn singular_seed_count(&self) -> usize {
        self.seeds
            .iter()
            .filter(|seed| matches!(seed.outcome(), TorsionXSeedLiftOutcome::SingularModuloPrime))
            .count()
    }

    #[cfg(test)]
    pub(super) fn uncertified_seed_count(&self) -> usize {
        self.seeds
            .iter()
            .filter(|seed| matches!(seed.outcome(), TorsionXSeedLiftOutcome::NotCertifiedInBound))
            .count()
    }
}

impl TorsionXSeedLiftReport {
    #[cfg(test)]
    pub(super) fn seed(&self) -> ReductionResidue {
        self.seed
    }

    pub(super) fn outcome(&self) -> &TorsionXSeedLiftOutcome {
        &self.outcome
    }
}

pub(super) fn lift_polynomial_roots_by_hensel(
    order: usize,
    source: TorsionXPolynomialSource,
    polynomial: &IntegerPolynomial,
    prime: ReductionPrime,
    seeds: Vec<ReductionResidue>,
) -> Result<TorsionXPolynomialLiftReport, HenselLiftError> {
    #[cfg(not(test))]
    let _ = (order, source);

    let root_bound = polynomial
        .cauchy_integer_root_bound()
        .ok_or(HenselLiftError::ConstantPolynomial)?;
    let prime_biguint = BigUint::from(prime.modulus());
    let seed_reports = seeds
        .into_iter()
        .map(|seed| {
            lift_seed_by_hensel(polynomial, prime, &prime_biguint, &root_bound, seed).map(
                |outcome| TorsionXSeedLiftReport {
                    #[cfg(test)]
                    seed,
                    outcome,
                },
            )
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(TorsionXPolynomialLiftReport {
        #[cfg(test)]
        order,
        #[cfg(test)]
        source,
        #[cfg(test)]
        prime,
        #[cfg(test)]
        root_bound,
        seeds: seed_reports,
    })
}

fn lift_seed_by_hensel(
    polynomial: &IntegerPolynomial,
    prime: ReductionPrime,
    prime_biguint: &BigUint,
    root_bound: &BigUint,
    seed: ReductionResidue,
) -> Result<TorsionXSeedLiftOutcome, HenselLiftError> {
    let seed_integer = BigInt::from(seed.representative());
    // check f′(x₀) ≡ 0 mod p
    let derivative_is_zero = prime
        .reduce_bigint(&polynomial.evaluate_derivative(&seed_integer))
        .is_zero();
    if derivative_is_zero {
        return Ok(TorsionXSeedLiftOutcome::SingularModuloPrime);
    }

    match hensel_lift_integer_root(polynomial, &seed_integer, prime_biguint, root_bound) {
        Ok(trace) => Ok(TorsionXSeedLiftOutcome::Certified(trace)),
        Err(HenselLiftError::IntegerRootNotCertifiedInBound) => {
            Ok(TorsionXSeedLiftOutcome::NotCertifiedInBound)
        }
        Err(error) => Err(error),
    }
}
