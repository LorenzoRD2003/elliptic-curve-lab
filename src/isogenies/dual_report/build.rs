use crate::isogenies::dual_report::summaries::{
    DegreeFactorizationSummary, KernelDescriptionSummary,
};
use num_bigint::BigUint;

/// One side of a duality report, either `phi` or `phi_hat`.
///
/// This groups together the ordinary degree, any finer
/// separable/inseparable factorization data, and the currently available
/// kernel summary for that same isogeny. Grouping them avoids parallel
/// `phi_*` / `dual_*` fields inside the main report.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DualIsogenySideSummary {
    degree: BigUint,
    degree_factorization: DegreeFactorizationSummary,
    kernel_summary: KernelDescriptionSummary,
}

impl DualIsogenySideSummary {
    pub(crate) fn new(
        degree: impl Into<BigUint>,
        degree_factorization: DegreeFactorizationSummary,
        kernel_summary: KernelDescriptionSummary,
    ) -> Self {
        Self {
            degree: degree.into(),
            degree_factorization,
            kernel_summary,
        }
    }

    pub(crate) fn degree(&self) -> &BigUint {
        &self.degree
    }

    pub(crate) fn degree_factorization(&self) -> &DegreeFactorizationSummary {
        &self.degree_factorization
    }

    pub(crate) fn kernel_summary(&self) -> &KernelDescriptionSummary {
        &self.kernel_summary
    }
}
