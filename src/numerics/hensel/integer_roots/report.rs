use num_bigint::{BigInt, BigUint};

use crate::numerics::hensel::{HenselIntegerRootSearchConfig, HenselLiftTrace};

/// Certified recovery of an integer root from a simple modular Hensel lift.
///
/// The report stores both sides of the certificate: the modular lift trace that
/// reached precision `pᵉ`, and the centered integer representative `r ∈ ℤ`
/// that passed the final exact check `f(r) = 0`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselIntegerRootTrace {
    lift_trace: HenselLiftTrace,
    root_bound: BigUint,
    modulus: BigUint,
    candidate_root: BigInt,
}

impl HenselIntegerRootTrace {
    pub(super) fn new(
        lift_trace: HenselLiftTrace,
        root_bound: BigUint,
        modulus: BigUint,
        candidate_root: BigInt,
    ) -> Self {
        Self {
            lift_trace,
            root_bound,
            modulus,
            candidate_root,
        }
    }

    /// Returns the underlying modular Hensel lift trace.
    pub(crate) fn lift_trace(&self) -> &HenselLiftTrace {
        &self.lift_trace
    }

    /// Returns the certified integer-root bound `B₀`.
    pub(crate) fn root_bound(&self) -> &BigUint {
        &self.root_bound
    }

    /// Returns the final modulus `pᵉ` used for unique recovery.
    pub(crate) fn modulus(&self) -> &BigUint {
        &self.modulus
    }

    /// Returns the certified integer root.
    pub(crate) fn candidate_root(&self) -> &BigInt {
        &self.candidate_root
    }
}

/// Report for one exhaustive simple-root seed scan modulo a chosen prime.
///
/// The counters distinguish the three outcomes that matter for this first
/// simple-root route:
///
/// 1. simple seeds that were lifted;
/// 2. singular modular roots that were intentionally skipped;
/// 3. simple seeds whose p-adic lift did not certify a bounded root in `ℤ`.
///
/// Complexity: the search routine sorts and deduplicates the `r` certified
/// roots before constructing this report, which costs `Θ(r log r)`. The report
/// itself stores `Θ(r)` roots and successful traces.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselIntegerRootSearchReport {
    config: HenselIntegerRootSearchConfig,
    simple_seed_count: usize,
    singular_seed_count: usize,
    uncertified_seed_count: usize,
    certified_roots: Vec<BigInt>,
    traces: Vec<HenselIntegerRootTrace>,
}

impl HenselIntegerRootSearchReport {
    pub(super) fn new(
        config: HenselIntegerRootSearchConfig,
        simple_seed_count: usize,
        singular_seed_count: usize,
        uncertified_seed_count: usize,
        certified_roots: Vec<BigInt>,
        traces: Vec<HenselIntegerRootTrace>,
    ) -> Self {
        Self {
            config,
            simple_seed_count,
            singular_seed_count,
            uncertified_seed_count,
            certified_roots,
            traces,
        }
    }

    /// Returns the config used for the scan.
    pub(crate) fn config(&self) -> &HenselIntegerRootSearchConfig {
        &self.config
    }

    /// Returns how many simple modular roots were lifted.
    pub(crate) fn simple_seed_count(&self) -> usize {
        self.simple_seed_count
    }

    /// Returns how many modular roots were skipped because the derivative
    /// vanished modulo `p`.
    pub(crate) fn singular_seed_count(&self) -> usize {
        self.singular_seed_count
    }

    /// Returns how many simple lifted seeds failed the final integer
    /// certification.
    pub(crate) fn uncertified_seed_count(&self) -> usize {
        self.uncertified_seed_count
    }

    /// Returns the certified bounded integer roots, sorted increasingly.
    pub(crate) fn certified_roots(&self) -> &[BigInt] {
        &self.certified_roots
    }

    /// Returns the successful lift traces, one per certified root.
    pub(crate) fn traces(&self) -> &[HenselIntegerRootTrace] {
        &self.traces
    }

    /// Returns whether the scan certified at least one integer root.
    pub(crate) fn has_certified_roots(&self) -> bool {
        !self.certified_roots.is_empty()
    }
}
