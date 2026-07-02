use num_bigint::BigUint;

/// One recorded discrete-log probe inside the `r`-Sylow subgroup.
///
/// The deliberately naive Exercise 3 route checks whether `α = xδ` for
/// `x = 1, …, r^k`. Each step stores the candidate `x` and the tested group
/// element `[x]δ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CyclicPrimeRootStep<P> {
    discrete_log_candidate: BigUint,
    candidate_multiple: P,
}

impl<P> CyclicPrimeRootStep<P> {
    pub(crate) fn new(discrete_log_candidate: BigUint, candidate_multiple: P) -> Self {
        Self {
            discrete_log_candidate,
            candidate_multiple,
        }
    }

    /// Returns the tested candidate `x`.
    pub(crate) fn discrete_log_candidate(&self) -> &BigUint {
        &self.discrete_log_candidate
    }

    /// Returns the tested group element `[x]δ`.
    pub(crate) fn candidate_multiple(&self) -> &P {
        &self.candidate_multiple
    }
}
