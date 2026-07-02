//! Perfect-power detection via bounded integer roots.

mod algorithm;
mod candidate;
mod exponents;
mod report;
mod validation;

#[cfg(test)]
mod tests;

pub(crate) use algorithm::detect_perfect_power;
pub(crate) use report::{PerfectPowerCandidateReport, PerfectPowerOutcome, PerfectPowerReport};

/// Configuration for the current Hensel-based perfect-power detector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PerfectPowerSearchConfig {
    max_seed_scan: u64,
}

impl Default for PerfectPowerSearchConfig {
    fn default() -> Self {
        Self {
            max_seed_scan: 10_000,
        }
    }
}

impl PerfectPowerSearchConfig {
    /// Sets the residue-scan limit passed to the Hensel integer-root search.
    ///
    /// Complexity: `Θ(1)`.
    pub(crate) fn with_max_seed_scan(mut self, max_seed_scan: u64) -> Self {
        self.max_seed_scan = max_seed_scan;
        self
    }

    /// Returns the residue-scan limit used for each modular root search.
    pub(crate) fn max_seed_scan(&self) -> u64 {
        self.max_seed_scan
    }
}
