use num_bigint::{BigInt, BigUint};

use crate::numerics::{hensel::HenselLiftError, perfect_powers::PerfectPowerSearchConfig};

/// Outcome of the Hensel-based perfect-power detector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PerfectPowerOutcome {
    /// A decomposition `N = aᵠ` was certified.
    PerfectPower {
        /// The certified base `a`.
        base: BigUint,
        /// The prime exponent `q`.
        exponent: u32,
    },
    /// No prime exponent `q ≤ ⌊log₂ N⌋` yielded a certified base.
    NotPerfectPower,
    /// The staged algorithm is only for `N > 1`.
    DegenerateInput,
    /// The staged algorithm assumes `gcd(N, 6) = 1`.
    NotCoprimeToSix,
    /// The input is too large for the current exponent/index representation.
    ExponentTooLarge,
    /// The underlying simple-root Hensel search could not run.
    HenselFailure(HenselLiftError),
}

/// Report for one prime-exponent candidate `q`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PerfectPowerCandidateReport {
    exponent: u32,
    hensel_prime: BigUint,
    root_bound: BigUint,
    certified_roots: Vec<BigInt>,
}

impl PerfectPowerCandidateReport {
    pub(super) fn new(
        exponent: u32,
        hensel_prime: BigUint,
        root_bound: BigUint,
        certified_roots: Vec<BigInt>,
    ) -> Self {
        Self {
            exponent,
            hensel_prime,
            root_bound,
            certified_roots,
        }
    }

    /// Returns the prime exponent `q` tested in `xᵠ − N`.
    pub(crate) fn exponent(&self) -> u32 {
        self.exponent
    }

    /// Returns the Hensel prime used for this exponent.
    pub(crate) fn hensel_prime(&self) -> &BigUint {
        &self.hensel_prime
    }

    /// Returns the bound used for candidate integer bases.
    pub(crate) fn root_bound(&self) -> &BigUint {
        &self.root_bound
    }

    /// Returns roots certified for `xᵠ − N` inside the bound.
    pub(crate) fn certified_roots(&self) -> &[BigInt] {
        &self.certified_roots
    }
}

/// Report for a perfect-power search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PerfectPowerReport {
    input: BigUint,
    config: PerfectPowerSearchConfig,
    candidate_reports: Vec<PerfectPowerCandidateReport>,
    outcome: PerfectPowerOutcome,
}

impl PerfectPowerReport {
    pub(super) fn new(
        input: BigUint,
        config: PerfectPowerSearchConfig,
        candidate_reports: Vec<PerfectPowerCandidateReport>,
        outcome: PerfectPowerOutcome,
    ) -> Self {
        Self {
            input,
            config,
            candidate_reports,
            outcome,
        }
    }

    /// Returns the tested input `N`.
    pub(crate) fn input(&self) -> &BigUint {
        &self.input
    }

    /// Returns the search configuration.
    pub(crate) fn config(&self) -> &PerfectPowerSearchConfig {
        &self.config
    }

    /// Returns the per-exponent candidate reports.
    pub(crate) fn candidate_reports(&self) -> &[PerfectPowerCandidateReport] {
        &self.candidate_reports
    }

    /// Returns the final outcome.
    pub(crate) fn outcome(&self) -> &PerfectPowerOutcome {
        &self.outcome
    }
}
