use num_complex::Complex64;

use crate::elliptic_curves::analytic::{
    AnalyticWeierstrassCurve,
    periods::{NumericalRecoveryMetadata, WeierstrassCubicRoots},
};
use crate::numerics::ComplexApproxComparison;

/// Structured report for one successful cubic-root recovery attempt.
#[derive(Clone, Debug, PartialEq)]
pub struct CubicRootRecoveryReport {
    curve: AnalyticWeierstrassCurve,
    roots: WeierstrassCubicRoots,
    g2_comparison: ComplexApproxComparison,
    g3_comparison: ComplexApproxComparison,
    metadata: NumericalRecoveryMetadata,
    cardano_diagnostics: Option<CardanoRootRecoveryDiagnostics>,
}

/// Root-recovery-specific Cardano-branch diagnostics.
///
/// These values describe the branch choice used during the current
/// Cardano/Newton hybrid route. They intentionally live on the cubic-root
/// recovery report rather than in the more general period-recovery metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct CardanoRootRecoveryDiagnostics {
    discriminant: Complex64,
    product_residual_norm: f64,
    selected_u_branch_index: usize,
    selected_v_branch_index: usize,
}

impl CubicRootRecoveryReport {
    pub(crate) fn new(
        curve: AnalyticWeierstrassCurve,
        roots: WeierstrassCubicRoots,
        g2_comparison: ComplexApproxComparison,
        g3_comparison: ComplexApproxComparison,
        metadata: NumericalRecoveryMetadata,
        cardano_diagnostics: Option<CardanoRootRecoveryDiagnostics>,
    ) -> Self {
        Self {
            curve,
            roots,
            g2_comparison,
            g3_comparison,
            metadata,
            cardano_diagnostics,
        }
    }

    /// Returns the original analytic curve.
    pub fn curve(&self) -> &AnalyticWeierstrassCurve {
        &self.curve
    }

    /// Returns the recovered cubic roots.
    pub fn roots(&self) -> &WeierstrassCubicRoots {
        &self.roots
    }

    /// Returns the comparison between reconstructed and original `g₂`.
    #[cfg(any(test, feature = "visualization"))]
    pub(crate) fn g2_comparison(&self) -> &ComplexApproxComparison {
        &self.g2_comparison
    }

    /// Returns the comparison between reconstructed and original `g₃`.
    #[cfg(any(test, feature = "visualization"))]
    pub(crate) fn g3_comparison(&self) -> &ComplexApproxComparison {
        &self.g3_comparison
    }

    /// Returns the numerical execution metadata for the recovery run.
    pub fn metadata(&self) -> &NumericalRecoveryMetadata {
        &self.metadata
    }

    /// Returns the Cardano-branch diagnostics, when the current recovery
    /// route recorded them.
    pub fn cardano_diagnostics(&self) -> Option<&CardanoRootRecoveryDiagnostics> {
        self.cardano_diagnostics.as_ref()
    }

    /// Returns the Cardano discriminant
    /// `(q/2)^2 + (p/3)^3`, when the run recorded it.
    pub fn cardano_discriminant(&self) -> Option<&Complex64> {
        self.cardano_diagnostics
            .as_ref()
            .map(CardanoRootRecoveryDiagnostics::discriminant)
    }

    /// Returns the residual norm of the selected Cardano branch condition
    /// `uv ≈ -p/3`, when the run recorded it.
    pub fn cardano_product_residual_norm(&self) -> Option<f64> {
        self.cardano_diagnostics
            .as_ref()
            .map(CardanoRootRecoveryDiagnostics::product_residual_norm)
    }

    /// Returns the selected branch index for `u`, when the run recorded it.
    pub fn selected_u_branch_index(&self) -> Option<usize> {
        self.cardano_diagnostics
            .as_ref()
            .map(CardanoRootRecoveryDiagnostics::selected_u_branch_index)
    }

    /// Returns the selected branch index for `v`, when the run recorded it.
    pub fn selected_v_branch_index(&self) -> Option<usize> {
        self.cardano_diagnostics
            .as_ref()
            .map(CardanoRootRecoveryDiagnostics::selected_v_branch_index)
    }

    /// Returns whether the selected Cardano pair used the principal
    /// cube-root branch for both `u` and `v`, when the run recorded the
    /// branch indices.
    pub fn used_principal_cardano_branches(&self) -> Option<bool> {
        self.cardano_diagnostics
            .as_ref()
            .map(CardanoRootRecoveryDiagnostics::used_principal_branches)
    }

    /// Returns the reconstructed `g₂` derived from the recovered roots.
    pub fn reconstructed_g2(&self) -> &Complex64 {
        self.g2_comparison.left()
    }

    /// Returns the original curve-side `g₂`.
    pub fn curve_g2(&self) -> &Complex64 {
        self.g2_comparison.right()
    }

    /// Returns the reconstructed `g₃` derived from the recovered roots.
    pub fn reconstructed_g3(&self) -> &Complex64 {
        self.g3_comparison.left()
    }

    /// Returns the original curve-side `g₃`.
    pub fn curve_g3(&self) -> &Complex64 {
        self.g3_comparison.right()
    }

    /// Returns whether both reconstructed coefficients agree approximately
    /// with the original curve-side coefficients.
    pub fn reconstruction_agrees(&self) -> bool {
        self.g2_comparison.agrees_approximately() && self.g3_comparison.agrees_approximately()
    }
}

impl CardanoRootRecoveryDiagnostics {
    pub(crate) fn new(
        discriminant: Complex64,
        product_residual_norm: f64,
        selected_u_branch_index: usize,
        selected_v_branch_index: usize,
    ) -> Self {
        Self {
            discriminant,
            product_residual_norm,
            selected_u_branch_index,
            selected_v_branch_index,
        }
    }

    /// Returns the Cardano discriminant `(q/2)^2 + (p/3)^3`.
    pub fn discriminant(&self) -> &Complex64 {
        &self.discriminant
    }

    /// Returns the residual norm of the selected branch condition
    /// `uv ≈ -p/3`.
    pub fn product_residual_norm(&self) -> f64 {
        self.product_residual_norm
    }

    /// Returns the selected branch index for `u`.
    pub fn selected_u_branch_index(&self) -> usize {
        self.selected_u_branch_index
    }

    /// Returns the selected branch index for `v`.
    pub fn selected_v_branch_index(&self) -> usize {
        self.selected_v_branch_index
    }

    /// Returns whether both selected indices are the principal branch.
    pub fn used_principal_branches(&self) -> bool {
        self.selected_u_branch_index == 0 && self.selected_v_branch_index == 0
    }
}
