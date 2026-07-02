use crate::elliptic_curves::analytic::{ModularMatrix, UpperHalfPlanePoint};

/// Records one actual modular transformation used while reducing `τ` toward
/// the standard fundamental domain.
#[derive(Clone, Debug, PartialEq)]
pub struct FundamentalDomainReductionStep {
    applied_matrix: ModularMatrix,
    before: UpperHalfPlanePoint,
    after: UpperHalfPlanePoint,
    reason: FundamentalDomainReductionStepReason,
}

impl FundamentalDomainReductionStep {
    pub(crate) fn new(
        applied_matrix: ModularMatrix,
        before: UpperHalfPlanePoint,
        after: UpperHalfPlanePoint,
        reason: FundamentalDomainReductionStepReason,
    ) -> Self {
        Self {
            applied_matrix,
            before,
            after,
            reason,
        }
    }

    /// Returns the matrix applied at this step.
    pub fn applied_matrix(&self) -> ModularMatrix {
        self.applied_matrix.clone()
    }

    /// Returns the upper-half-plane point before the step.
    pub fn before(&self) -> &UpperHalfPlanePoint {
        &self.before
    }

    /// Returns the upper-half-plane point after the step.
    pub fn after(&self) -> &UpperHalfPlanePoint {
        &self.after
    }

    /// Returns why this modular transformation was applied.
    pub fn reason(&self) -> FundamentalDomainReductionStepReason {
        self.reason
    }
}

/// Explains why one modular step was applied during reduction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FundamentalDomainReductionStepReason {
    /// A translation `τ ↦ τ - n` was needed because `Re(τ)` lay outside the
    /// centered strip `|Re(τ)| ≤ 1/2`.
    RealPartOutsideCenteredStrip,
    /// An inversion `τ ↦ -1/τ` was needed because `|τ| < 1`.
    NormLessThanOne,
}

/// Final status of one attempted reduction to the standard fundamental domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FundamentalDomainReductionStatus {
    /// The input already satisfied the standard domain conditions, so no
    /// modular step was needed.
    AlreadyReduced,
    /// The algorithm reached a point satisfying the current domain test.
    Reduced,
    /// The algorithm stopped after `max_steps` transformations without yet
    /// reaching the target region.
    StepLimitReached,
}

/// Structured report for reducing `τ` to the classical standard fundamental
/// domain of `SL_2(ℤ)`.
#[derive(Clone, Debug, PartialEq)]
pub struct FundamentalDomainReductionReport {
    original_tau: UpperHalfPlanePoint,
    reduced_tau: UpperHalfPlanePoint,
    accumulated_matrix: ModularMatrix,
    steps: Vec<FundamentalDomainReductionStep>,
    status: FundamentalDomainReductionStatus,
}

impl FundamentalDomainReductionReport {
    pub(crate) fn new(
        original_tau: UpperHalfPlanePoint,
        reduced_tau: UpperHalfPlanePoint,
        accumulated_matrix: ModularMatrix,
        steps: Vec<FundamentalDomainReductionStep>,
        status: FundamentalDomainReductionStatus,
    ) -> Self {
        Self {
            original_tau,
            reduced_tau,
            accumulated_matrix,
            steps,
            status,
        }
    }

    /// Returns the original upper-half-plane input.
    pub fn original_tau(&self) -> &UpperHalfPlanePoint {
        &self.original_tau
    }

    /// Returns the final representative reached by the reduction routine.
    pub fn reduced_tau(&self) -> &UpperHalfPlanePoint {
        &self.reduced_tau
    }

    /// Returns the accumulated modular matrix `γ` such that the final point is
    /// `γ(τ)` for the original input `τ`.
    pub fn accumulated_matrix(&self) -> ModularMatrix {
        self.accumulated_matrix.clone()
    }

    /// Returns the sequence of actual modular transformations that were
    /// applied.
    pub fn steps(&self) -> &[FundamentalDomainReductionStep] {
        &self.steps
    }

    /// Returns the terminal status of the reduction attempt.
    pub fn status(&self) -> FundamentalDomainReductionStatus {
        self.status
    }

    /// Returns whether the final point was accepted as lying in the target
    /// standard fundamental domain.
    pub fn is_reduced(&self) -> bool {
        matches!(
            self.status,
            FundamentalDomainReductionStatus::AlreadyReduced
                | FundamentalDomainReductionStatus::Reduced
        )
    }
}
