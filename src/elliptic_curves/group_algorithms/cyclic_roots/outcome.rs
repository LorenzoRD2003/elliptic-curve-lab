/// Outcome of one attempted equation `[r]ρ = γ` in a cyclic group.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CyclicPrimeRootOutcome<P> {
    /// A root `ρ` was found.
    Root { root: P },
    /// The input `γ` has no `r`-th root in the represented cyclic group.
    NoRoot,
}

impl<P> CyclicPrimeRootOutcome<P> {
    /// Returns the root `ρ`, if one was found.
    pub fn root(&self) -> Option<&P> {
        match self {
            Self::Root { root } => Some(root),
            Self::NoRoot => None,
        }
    }

    /// Returns whether the attempt found an `r`-th root.
    pub fn found_root(&self) -> bool {
        matches!(self, Self::Root { .. })
    }
}
