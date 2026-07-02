use num_bigint::BigUint;

use crate::elliptic_curves::group_algorithms::cyclic_roots::{
    CyclicPrimeRootBezout, CyclicPrimeRootStep,
};

/// Route trace for one attempted prime-degree root extraction.
///
/// The trace records the group elements named in the exercise:
///
/// - `α = aγ`
/// - `β = r^k γ`
/// - a brute-force search path for `x` with `α = xδ`
/// - optional Bezout data for `s a + t r^(k+1) = 1`
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CyclicPrimeRootTrace<P> {
    alpha: Option<P>,
    beta: Option<P>,
    discrete_log: Option<BigUint>,
    bezout: Option<CyclicPrimeRootBezout>,
    steps: Vec<CyclicPrimeRootStep<P>>,
}

impl<P> CyclicPrimeRootTrace<P> {
    pub(crate) fn new(
        alpha: Option<P>,
        beta: Option<P>,
        discrete_log: Option<BigUint>,
        bezout: Option<CyclicPrimeRootBezout>,
        steps: Vec<CyclicPrimeRootStep<P>>,
    ) -> Self {
        Self {
            alpha,
            beta,
            discrete_log,
            bezout,
            steps,
        }
    }

    /// Returns `α = aγ`, once computed.
    pub fn alpha(&self) -> Option<&P> {
        self.alpha.as_ref()
    }

    /// Returns `β = r^k γ`, once computed.
    pub fn beta(&self) -> Option<&P> {
        self.beta.as_ref()
    }

    /// Returns the discovered discrete logarithm `x` with `α = xδ`.
    pub fn discrete_log(&self) -> Option<&BigUint> {
        self.discrete_log.as_ref()
    }

    /// Returns the Bezout data used in the final root formula.
    pub fn bezout(&self) -> Option<&CyclicPrimeRootBezout> {
        self.bezout.as_ref()
    }

    /// Returns the brute-force discrete-log probes in order.
    pub fn steps(&self) -> &[CyclicPrimeRootStep<P>] {
        &self.steps
    }
}
