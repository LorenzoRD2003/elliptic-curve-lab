use crate::elliptic_curves::frobenius::HasseInterval;

/// One tested candidate `M ∈ H(q)` in the naive Hasse-interval search.
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

/// Report for the naive search of an annihilating multiple inside `H(q)`.
///
/// Starting from the lower endpoint of the discrete Hasse interval, this
/// report records the tested sequence
///
/// `[L]P, [(L+1)]P, ..., [U]P`
///
/// where `H(q) = [L, U]`, until the first candidate `M` with `[M]P = O` is
/// found or the interval is exhausted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HasseMultipleSearchReport<P> {
    q: u128,
    interval: HasseInterval,
    tested_candidates: u128,
    first_annihilating_multiple: Option<u128>,
    steps: Vec<HasseMultipleSearchStep<P>>,
}

impl<P> HasseMultipleSearchReport<P> {
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

pub(crate) fn hasse_multiple_search_report<P>(
    q: u128,
    interval: HasseInterval,
    first_annihilating_multiple: Option<u128>,
    steps: Vec<HasseMultipleSearchStep<P>>,
) -> HasseMultipleSearchReport<P> {
    HasseMultipleSearchReport {
        q,
        interval,
        tested_candidates: steps.len() as u128,
        first_annihilating_multiple,
        steps,
    }
}
