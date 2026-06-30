use num_bigint::{BigInt, BigUint};

/// One successful simple-root Hensel step `x_{k+1} = x_k + t p^k`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselLiftStep {
    level: u32,
    root_before: BigInt,
    lift_digit: BigUint,
    root_after: BigInt,
}

impl HenselLiftStep {
    pub(super) fn new(
        level: u32,
        root_before: BigInt,
        lift_digit: BigUint,
        root_after: BigInt,
    ) -> Self {
        Self {
            level,
            root_before,
            lift_digit,
            root_after,
        }
    }

    /// Returns the source precision `k` in the step from modulo `p^k` to
    /// modulo `p^(k+1)`.
    pub(crate) fn level(&self) -> u32 {
        self.level
    }

    /// Returns the root representative before the step.
    pub(crate) fn root_before(&self) -> &BigInt {
        &self.root_before
    }

    /// Returns the digit `t` chosen modulo `p`.
    pub(crate) fn lift_digit(&self) -> &BigUint {
        &self.lift_digit
    }

    /// Returns the root representative after the step.
    pub(crate) fn root_after(&self) -> &BigInt {
        &self.root_after
    }
}

/// One successful fast square-root Hensel step from `p^k` to `p^m`, where
/// `m = min(2k, target_level)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselSquareRootFastStep {
    source_level: u32,
    target_level: u32,
    root_before: BigInt,
    root_after: BigInt,
}

impl HenselSquareRootFastStep {
    pub(super) fn new(
        source_level: u32,
        target_level: u32,
        root_before: BigInt,
        root_after: BigInt,
    ) -> Self {
        Self {
            source_level,
            target_level,
            root_before,
            root_after,
        }
    }

    /// Returns the source precision `k`.
    pub(crate) fn source_level(&self) -> u32 {
        self.source_level
    }

    /// Returns the target precision reached by this Newton-Hensel step.
    pub(crate) fn target_level(&self) -> u32 {
        self.target_level
    }

    /// Returns the root representative before the step.
    pub(crate) fn root_before(&self) -> &BigInt {
        &self.root_before
    }

    /// Returns the root representative after the step.
    pub(crate) fn root_after(&self) -> &BigInt {
        &self.root_after
    }
}
