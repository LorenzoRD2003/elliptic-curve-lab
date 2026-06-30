use num_bigint::{BigInt, BigUint};

use super::step::{HenselLiftStep, HenselSquareRootFastStep};

/// Educational trace for repeated simple-root Hensel lifting.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselLiftTrace {
    prime: BigUint,
    coefficients: Vec<BigInt>,
    initial_root: BigInt,
    steps: Vec<HenselLiftStep>,
}

impl HenselLiftTrace {
    pub(super) fn new(
        prime: BigUint,
        coefficients: Vec<BigInt>,
        initial_root: BigInt,
        steps: Vec<HenselLiftStep>,
    ) -> Self {
        Self {
            prime,
            coefficients,
            initial_root,
            steps,
        }
    }

    /// Returns the prime `p`.
    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the lifted polynomial coefficients in ascending degree order.
    pub(crate) fn coefficients(&self) -> &[BigInt] {
        &self.coefficients
    }

    /// Returns the initial root representative modulo `p`.
    pub(crate) fn initial_root(&self) -> &BigInt {
        &self.initial_root
    }

    /// Returns the recorded Hensel steps.
    pub(crate) fn steps(&self) -> &[HenselLiftStep] {
        &self.steps
    }

    /// Returns the final lifted root representative.
    pub(crate) fn final_root(&self) -> &BigInt {
        self.steps
            .last()
            .map_or(&self.initial_root, HenselLiftStep::root_after)
    }

    /// Returns the precision level reached by the trace.
    pub(crate) fn reached_level(&self) -> u32 {
        self.steps.last().map_or(1, |step| step.level() + 1)
    }
}

/// Educational trace for fast square-root Hensel lifting.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct HenselSquareRootFastTrace {
    prime: BigUint,
    value: BigInt,
    initial_root: BigInt,
    target_level: u32,
    steps: Vec<HenselSquareRootFastStep>,
}

impl HenselSquareRootFastTrace {
    pub(super) fn new(
        prime: BigUint,
        value: BigInt,
        initial_root: BigInt,
        target_level: u32,
        steps: Vec<HenselSquareRootFastStep>,
    ) -> Self {
        Self {
            prime,
            value,
            initial_root,
            target_level,
            steps,
        }
    }

    /// Returns the prime `p`.
    pub(crate) fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the radicand `a` in `x^2 = a`.
    pub(crate) fn value(&self) -> &BigInt {
        &self.value
    }

    /// Returns the initial root representative modulo `p`.
    pub(crate) fn initial_root(&self) -> &BigInt {
        &self.initial_root
    }

    /// Returns the requested final precision level.
    pub(crate) fn target_level(&self) -> u32 {
        self.target_level
    }

    /// Returns the recorded fast Newton-Hensel steps.
    pub(crate) fn steps(&self) -> &[HenselSquareRootFastStep] {
        &self.steps
    }

    /// Returns the final lifted root representative.
    pub(crate) fn final_root(&self) -> &BigInt {
        self.steps
            .last()
            .map_or(&self.initial_root, HenselSquareRootFastStep::root_after)
    }

    /// Returns the precision level reached by the trace.
    pub(crate) fn reached_level(&self) -> u32 {
        self.steps
            .last()
            .map_or(1, HenselSquareRootFastStep::target_level)
    }
}
