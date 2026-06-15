use crate::isogenies::dual_report::summaries::{
    DegreeFactorizationSummary, KernelDescriptionSummary,
};

/// One side of a duality report, either `phi` or `phi_hat`.
///
/// This groups together the ordinary degree, any finer
/// separable/inseparable factorization data, and the currently available
/// kernel summary for that same isogeny. Grouping them avoids parallel
/// `phi_*` / `dual_*` fields inside the main report.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DualIsogenySideSummary {
    degree: usize,
    degree_factorization: DegreeFactorizationSummary,
    kernel_summary: KernelDescriptionSummary,
}

impl DualIsogenySideSummary {
    pub(crate) fn new(
        degree: usize,
        degree_factorization: DegreeFactorizationSummary,
        kernel_summary: KernelDescriptionSummary,
    ) -> Self {
        Self {
            degree,
            degree_factorization,
            kernel_summary,
        }
    }

    pub(crate) fn degree(&self) -> usize {
        self.degree
    }

    pub(crate) fn degree_factorization(&self) -> &DegreeFactorizationSummary {
        &self.degree_factorization
    }

    pub(crate) fn kernel_summary(&self) -> &KernelDescriptionSummary {
        &self.kernel_summary
    }
}
