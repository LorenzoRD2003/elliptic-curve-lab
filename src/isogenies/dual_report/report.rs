use core::marker::PhantomData;

use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::dual_report::{
    build::DualIsogenySideSummary,
    summaries::{DegreeFactorizationSummary, DualityKind, KernelDescriptionSummary},
};

/// Structured duality report for two related isogeny objects.
///
/// This report is intentionally lightweight. It does not duplicate full curve
/// or map data. Instead, it records:
///
/// - the kind of duality story currently being modeled
/// - one summary for `phi`
/// - one summary for `phi_hat`
/// - whether the left and right duality relations currently verify
/// - free-form notes explaining the current mathematical scope
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DualIsogenyReport<Domain: CurveModel, Codomain: CurveModel> {
    duality_kind: DualityKind,
    phi: DualIsogenySideSummary,
    dual: DualIsogenySideSummary,
    left_relation_holds: bool,
    right_relation_holds: bool,
    notes: Vec<String>,
    marker: PhantomData<(Domain, Codomain)>,
}

impl<Domain: CurveModel, Codomain: CurveModel> DualIsogenyReport<Domain, Codomain> {
    pub(crate) fn new(
        duality_kind: DualityKind,
        phi: DualIsogenySideSummary,
        dual: DualIsogenySideSummary,
        left_relation_holds: bool,
        right_relation_holds: bool,
        notes: Vec<String>,
    ) -> Self {
        Self {
            duality_kind,
            phi,
            dual,
            left_relation_holds,
            right_relation_holds,
            notes,
            marker: PhantomData,
        }
    }

    /// Returns the high-level duality kind modeled by this report.
    pub fn duality_kind(&self) -> DualityKind {
        self.duality_kind
    }

    /// Returns `deg(phi)`.
    pub fn phi_degree(&self) -> usize {
        self.phi.degree()
    }

    /// Returns `deg(phi_hat)`.
    pub fn dual_degree(&self) -> usize {
        self.dual.degree()
    }

    /// Returns the separable/inseparable degree summary for `phi`.
    pub fn phi_degree_factorization(&self) -> &DegreeFactorizationSummary {
        self.phi.degree_factorization()
    }

    /// Returns the separable/inseparable degree summary for `phi_hat`.
    pub fn dual_degree_factorization(&self) -> &DegreeFactorizationSummary {
        self.dual.degree_factorization()
    }

    /// Returns the kernel summary for `phi`.
    pub fn phi_kernel_summary(&self) -> &KernelDescriptionSummary {
        self.phi.kernel_summary()
    }

    /// Returns the kernel summary for `phi_hat`.
    pub fn dual_kernel_summary(&self) -> &KernelDescriptionSummary {
        self.dual.kernel_summary()
    }

    /// Returns whether the left duality relation currently verifies.
    pub fn left_relation_holds(&self) -> bool {
        self.left_relation_holds
    }

    /// Returns whether the right duality relation currently verifies.
    pub fn right_relation_holds(&self) -> bool {
        self.right_relation_holds
    }

    /// Returns explanatory notes about the current mathematical scope.
    pub fn notes(&self) -> &[String] {
        &self.notes
    }
}
