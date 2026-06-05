use super::super::AnalyticCurveError;
use crate::numerics::ApproxTolerance;

/// Validated numerical policy bundle for period recovery and inverse
/// uniformization experiments.
///
/// The current knobs are intentionally small and explicit:
/// - one shared comparison tolerance
/// - iteration budgets for Newton and AGM phases
/// - a step budget for Abel-Jacobi or elliptic-integral quadrature
/// - a finite search radius over nearby lattice branches
///
/// All integer budgets must be strictly positive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeriodRecoveryConfig {
    tolerance: ApproxTolerance,
    newton_max_iterations: usize,
    agm_max_iterations: usize,
    abel_jacobi_integration_steps: usize,
    branch_lattice_search_radius: usize,
}

impl PeriodRecoveryConfig {
    /// Builds a validated period-recovery configuration from explicit knobs.
    pub fn new(
        tolerance: ApproxTolerance,
        newton_max_iterations: usize,
        agm_max_iterations: usize,
        abel_jacobi_integration_steps: usize,
        branch_lattice_search_radius: usize,
    ) -> Result<Self, AnalyticCurveError> {
        if newton_max_iterations == 0
            || agm_max_iterations == 0
            || abel_jacobi_integration_steps == 0
            || branch_lattice_search_radius == 0
        {
            return Err(AnalyticCurveError::InvalidPeriodRecoveryConfig);
        }

        Ok(Self {
            tolerance,
            newton_max_iterations,
            agm_max_iterations,
            abel_jacobi_integration_steps,
            branch_lattice_search_radius,
        })
    }

    /// Returns the tolerance policy used for numerical comparisons.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the Newton-iteration budget.
    pub fn newton_max_iterations(&self) -> usize {
        self.newton_max_iterations
    }

    /// Returns the AGM-iteration budget.
    pub fn agm_max_iterations(&self) -> usize {
        self.agm_max_iterations
    }

    /// Returns the quadrature step budget for Abel-Jacobi style integration.
    pub fn abel_jacobi_integration_steps(&self) -> usize {
        self.abel_jacobi_integration_steps
    }

    /// Returns the search radius over nearby lattice branches.
    pub fn branch_lattice_search_radius(&self) -> usize {
        self.branch_lattice_search_radius
    }

    /// Returns the baseline preset for educational experiments.
    pub fn educational_default() -> Self {
        Self::new(ApproxTolerance::educational_default(), 12, 10, 256, 2)
            .expect("educational period-recovery preset must stay valid")
    }

    /// Returns a tighter preset for more delicate recovery experiments.
    pub fn strict() -> Self {
        Self::new(ApproxTolerance::strict(), 20, 16, 512, 4)
            .expect("strict period-recovery preset must stay valid")
    }

    /// Returns a more permissive preset for coarse exploratory work.
    pub fn loose() -> Self {
        Self::new(ApproxTolerance::loose(), 8, 6, 128, 1)
            .expect("loose period-recovery preset must stay valid")
    }
}
