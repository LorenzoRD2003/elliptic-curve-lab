use crate::elliptic_curves::frobenius::HasseInterval;

/// One tested candidate `M ∈ H(q)` in one Hasse-interval multiple search.
///
/// The stored image is the actual group element `[M]P`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseMultipleSearchStep<P> {
    candidate_multiple: u128,
    image: P,
}

impl<P> HasseMultipleSearchStep<P> {
    pub(crate) fn new(candidate_multiple: u128, image: P) -> Self {
        Self {
            candidate_multiple,
            image,
        }
    }

    /// Returns the tested candidate `M`.
    pub fn candidate_multiple(&self) -> u128 {
        self.candidate_multiple
    }

    /// Returns the group element `[M]P`.
    pub fn image(&self) -> &P {
        &self.image
    }
}

/// Report for one search of an annihilating multiple inside `H(q)`.
///
/// The report records the tested sequence of integer candidates and the
/// corresponding group elements `[M]P` in the order the chosen search engine
/// explored them.
///
/// For the current naive route that order is left-to-right from the lower
/// endpoint. Future or alternative search engines may explore the same
/// interval in a different order while still reporting their tested path
/// through this same type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseMultipleSearchReport<P> {
    q: u128,
    interval: HasseInterval,
    tested_candidates: u128,
    first_annihilating_multiple: Option<u128>,
    steps: Vec<HasseMultipleSearchStep<P>>,
}

impl<P> HasseMultipleSearchReport<P> {
    pub(crate) fn new(
        q: u128,
        interval: HasseInterval,
        first_annihilating_multiple: Option<u128>,
        steps: Vec<HasseMultipleSearchStep<P>>,
    ) -> Self {
        Self {
            q,
            interval,
            tested_candidates: steps.len() as u128,
            first_annihilating_multiple,
            steps,
        }
    }

    /// Returns the finite field order `q`.
    pub fn q(&self) -> u128 {
        self.q
    }

    /// Returns the searched Hasse interval `H(q)`.
    pub fn interval(&self) -> &HasseInterval {
        &self.interval
    }

    /// Returns how many interval candidates were tested.
    pub fn tested_candidates(&self) -> u128 {
        self.tested_candidates
    }

    /// Returns the first `M ∈ H(q)` such that `[M]P = O`, if found.
    pub fn first_annihilating_multiple(&self) -> Option<u128> {
        self.first_annihilating_multiple
    }

    /// Returns whether the search found an annihilating multiple.
    pub fn found(&self) -> bool {
        self.first_annihilating_multiple.is_some()
    }

    /// Returns the tested steps in order.
    pub fn steps(&self) -> &[HasseMultipleSearchStep<P>] {
        &self.steps
    }
}
