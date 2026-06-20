use crate::elliptic_curves::frobenius::schoof::{
    SchoofTraceMod2Report, SchoofTraceModOddPrimeReport,
};
use crate::fields::traits::FiniteField;
use crate::numerics::chinese_remainder::{ChineseRemainderSolution, Congruence};

/// Final outcome of the current CRT recombination stage in Schoof's algorithm.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SchoofTraceCrtOutcome {
    /// Every requested prime step contributed one trace congruence, so CRT
    /// produced one combined class modulo the product of those moduli.
    Combined { solution: ChineseRemainderSolution },
    /// One odd-prime step stopped before producing a congruence, so the CRT
    /// stage returns the partial solution accumulated strictly before that
    /// prime.
    BlockedOnOddPrime {
        blocked_prime: usize,
        partial_solution: ChineseRemainderSolution,
    },
}

/// Report for CRT recombination across the currently available Schoof prime steps.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofTraceCrtReport<F: FiniteField> {
    field_order: u128,
    mod_2_report: SchoofTraceMod2Report<F>,
    odd_prime_reports: Vec<SchoofTraceModOddPrimeReport<F>>,
    resolved_congruences: Vec<Congruence>,
    outcome: SchoofTraceCrtOutcome,
}

impl<F: FiniteField> SchoofTraceCrtReport<F> {
    pub(crate) fn new(
        field_order: u128,
        mod_2_report: SchoofTraceMod2Report<F>,
        odd_prime_reports: Vec<SchoofTraceModOddPrimeReport<F>>,
        resolved_congruences: Vec<Congruence>,
        outcome: SchoofTraceCrtOutcome,
    ) -> Self {
        Self {
            field_order,
            mod_2_report,
            odd_prime_reports,
            resolved_congruences,
            outcome,
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        self.field_order
    }

    /// Returns the initial `ℓ = 2` Schoof report.
    pub fn mod_2_report(&self) -> &SchoofTraceMod2Report<F> {
        &self.mod_2_report
    }

    /// Returns the odd-prime step reports in the order they were attempted.
    pub fn odd_prime_reports(&self) -> &[SchoofTraceModOddPrimeReport<F>] {
        &self.odd_prime_reports
    }

    /// Returns the congruences that successfully entered the CRT stage.
    pub fn resolved_congruences(&self) -> &[Congruence] {
        &self.resolved_congruences
    }

    /// Returns the current CRT-stage outcome.
    pub fn outcome(&self) -> &SchoofTraceCrtOutcome {
        &self.outcome
    }

    /// Returns the combined CRT class when every requested prime step resolved.
    pub fn combined_solution(&self) -> Option<&ChineseRemainderSolution> {
        match &self.outcome {
            SchoofTraceCrtOutcome::Combined { solution } => Some(solution),
            SchoofTraceCrtOutcome::BlockedOnOddPrime { .. } => None,
        }
    }
}
