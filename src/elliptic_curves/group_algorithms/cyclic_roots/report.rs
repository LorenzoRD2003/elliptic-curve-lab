use crate::elliptic_curves::group_algorithms::cyclic_roots::{
    CyclicPrimeRootInput, CyclicPrimeRootOutcome, CyclicPrimeRootTrace,
};

/// Full report for one prime-degree root extraction attempt.
///
/// This is intentionally route-preserving: future implementations should store
/// the original target `γ`, the supplied `r`-Sylow generator `δ`, the integer
/// decomposition `|G| = a r^k`, and the intermediate group elements that lead
/// to the final outcome.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CyclicPrimeRootReport<P> {
    input: CyclicPrimeRootInput,
    target: P,
    sylow_generator: P,
    trace: CyclicPrimeRootTrace<P>,
    outcome: CyclicPrimeRootOutcome<P>,
}

impl<P> CyclicPrimeRootReport<P> {
    pub(crate) fn new(
        input: CyclicPrimeRootInput,
        target: P,
        sylow_generator: P,
        trace: CyclicPrimeRootTrace<P>,
        outcome: CyclicPrimeRootOutcome<P>,
    ) -> Self {
        Self {
            input,
            target,
            sylow_generator,
            trace,
            outcome,
        }
    }

    /// Returns the integer-side setup for this root attempt.
    pub(crate) fn input(&self) -> &CyclicPrimeRootInput {
        &self.input
    }

    /// Returns the target group element `γ`.
    pub(crate) fn target(&self) -> &P {
        &self.target
    }

    /// Returns the supplied `r`-Sylow generator `δ`.
    pub(crate) fn sylow_generator(&self) -> &P {
        &self.sylow_generator
    }

    /// Returns the route trace.
    pub(crate) fn trace(&self) -> &CyclicPrimeRootTrace<P> {
        &self.trace
    }

    /// Returns the final outcome.
    pub(crate) fn outcome(&self) -> &CyclicPrimeRootOutcome<P> {
        &self.outcome
    }

    /// Returns the root `ρ`, if one was found.
    pub(crate) fn root(&self) -> Option<&P> {
        self.outcome.root()
    }
}
