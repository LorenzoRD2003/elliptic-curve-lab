use num_bigint::BigUint;

/// One prime-by-prime reduction step in the order-from-multiple algorithm.
///
/// If the supplied multiple has the factorization `M = Π ℓᵢ^eᵢ`, the
/// short-Weierstrass wrapper isolates one prime-primary component at a time
/// and then records how much of that `ℓ`-power can be removed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointOrderReductionStep {
    prime: BigUint,
    exponent_in_multiple: u32,
    removed_exponent: u32,
    remaining_multiple_after_step: BigUint,
}

impl PointOrderReductionStep {
    pub fn new(
        prime: BigUint,
        exponent_in_multiple: u32,
        removed_exponent: u32,
        remaining_multiple_after_step: BigUint,
    ) -> Self {
        Self {
            prime,
            exponent_in_multiple,
            removed_exponent,
            remaining_multiple_after_step,
        }
    }

    /// Returns the prime `ℓ` considered at this step.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the exponent `e` of `ℓ^e` in the supplied multiple.
    pub fn exponent_in_multiple(&self) -> u32 {
        self.exponent_in_multiple
    }

    /// Returns how many copies of `ℓ` were removed while preserving
    /// annihilation of the point.
    pub fn removed_exponent(&self) -> u32 {
        self.removed_exponent
    }

    /// Returns the remaining local exponent after this prime reduction.
    pub fn remaining_exponent(&self) -> u32 {
        self.exponent_in_multiple - self.removed_exponent
    }

    /// Returns the running remaining multiple after finishing this prime.
    pub fn remaining_multiple_after_step(&self) -> &BigUint {
        &self.remaining_multiple_after_step
    }
}

/// Report for recovering the exact order of a point from one known multiple.
///
/// Starting from one annihilating multiple `M` with `[M]P = O`, factor
/// `M = Π ℓᵢ^eᵢ` and isolate each `ℓ`-primary component. The final
/// reconstructed product of the local exact powers is the exact order of `P`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointOrderFromMultipleReport {
    supplied_multiple: BigUint,
    exact_order: BigUint,
    steps: Vec<PointOrderReductionStep>,
}

impl PointOrderFromMultipleReport {
    pub fn new(
        supplied_multiple: BigUint,
        exact_order: BigUint,
        steps: Vec<PointOrderReductionStep>,
    ) -> Self {
        Self {
            supplied_multiple,
            exact_order,
            steps,
        }
    }

    /// Returns the original supplied multiple `M`.
    pub fn supplied_multiple(&self) -> &BigUint {
        &self.supplied_multiple
    }

    /// Returns the recovered exact order of the point.
    pub fn exact_order(&self) -> &BigUint {
        &self.exact_order
    }

    /// Returns the remaining multiple after all prime reductions.
    ///
    /// In the current algorithm this equals [`Self::exact_order`]. It remains
    /// as a derived accessor because some educational visualizations phrase the
    /// final answer as “the remaining multiple after all reductions”.
    pub fn remaining_multiple(&self) -> &BigUint {
        self.exact_order()
    }

    /// Returns the per-prime reduction steps.
    pub fn steps(&self) -> &[PointOrderReductionStep] {
        &self.steps
    }
}
