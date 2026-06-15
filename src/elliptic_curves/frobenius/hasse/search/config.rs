/// Traversal policy for the Hasse-interval baby-step/giant-step search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum HasseBsgsTraversal {
    /// The current optimized default.
    LeftToRight,
    /// Center-first traversal inspired by Lecture 7, §7.11.
    ///
    /// The heuristic motivation is distributional: if the decisive
    /// annihilating multiple `M` is statistically more likely to lie near the
    /// center `q + 1` of `H(q)`, then a traversal that probes central blocks
    /// first may terminate earlier on average than a monotone left-to-right
    /// sweep.
    ///
    /// The current implementation avoids the earlier long-jump pathology by
    /// keeping two giant-step frontiers, one expanding to the right of the
    /// center block and one to the left, so each newly visited block costs
    /// only one giant-step update.
    MiddleOut,
}

/// Optional parity information for the unknown annihilating multiple.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum HasseBsgsParity {
    Unknown,
    Even,
    Odd,
}

/// Internal configuration for Hasse-interval BSGS.
///
/// This type is intentionally small and builder-like:
///
/// - [`Default`] / [`Self::new`] define the crate's current tuned baseline
/// - `with_...` methods let benchmarks and internal experiments vary one knob
///   at a time without exposing struct literals at call sites
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct HasseBsgsConfig {
    traversal: HasseBsgsTraversal,
    use_fast_negation: bool,
    known_parity: HasseBsgsParity,
}

impl HasseBsgsConfig {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    #[allow(dead_code)]
    pub(crate) fn with_traversal(mut self, traversal: HasseBsgsTraversal) -> Self {
        self.traversal = traversal;
        self
    }

    #[allow(dead_code)]
    pub(crate) fn with_fast_negation(mut self, use_fast_negation: bool) -> Self {
        self.use_fast_negation = use_fast_negation;
        self
    }

    pub(crate) fn with_known_parity(mut self, known_parity: HasseBsgsParity) -> Self {
        self.known_parity = known_parity;
        self
    }

    pub(crate) fn traversal(self) -> HasseBsgsTraversal {
        self.traversal
    }

    pub(crate) fn uses_fast_negation(self) -> bool {
        self.use_fast_negation
    }

    pub(crate) fn known_parity(self) -> HasseBsgsParity {
        self.known_parity
    }
}

impl Default for HasseBsgsConfig {
    fn default() -> Self {
        Self {
            traversal: HasseBsgsTraversal::LeftToRight,
            use_fast_negation: true,
            known_parity: HasseBsgsParity::Unknown,
        }
    }
}
