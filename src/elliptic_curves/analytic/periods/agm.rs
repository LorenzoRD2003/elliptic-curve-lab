use num_complex::Complex64;

use crate::elliptic_curves::analytic::AnalyticCurveError;
use crate::elliptic_curves::analytic::periods::config::PeriodRecoveryConfig;
use crate::fields::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Validated numerical policy for the raw complex AGM iteration.
///
/// This value object is intentionally narrower than [`PeriodRecoveryConfig`]:
/// the raw AGM primitive only needs one comparison tolerance and one iteration
/// budget. Higher-level period-recovery or elliptic-integral layers can derive
/// this smaller config from a broader recovery policy when needed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexAgmConfig {
    tolerance: ApproxTolerance,
    max_iterations: usize,
}

impl ComplexAgmConfig {
    /// Builds a validated complex-AGM policy from explicit knobs.
    ///
    /// The tolerance components must be finite, nonnegative, and not both
    /// zero. The iteration budget must be strictly positive.
    pub fn new(
        tolerance: ApproxTolerance,
        max_iterations: usize,
    ) -> Result<Self, AnalyticCurveError> {
        let tolerance_is_valid = tolerance.absolute.is_finite()
            && tolerance.relative.is_finite()
            && tolerance.absolute >= 0.0
            && tolerance.relative >= 0.0
            && (tolerance.absolute > 0.0 || tolerance.relative > 0.0);

        if !tolerance_is_valid || max_iterations == 0 {
            return Err(AnalyticCurveError::InvalidPeriodRecoveryConfig);
        }

        Ok(Self {
            tolerance,
            max_iterations,
        })
    }

    /// Derives the raw AGM policy from a broader period-recovery config.
    ///
    /// This keeps only the tolerance and AGM iteration budget, intentionally
    /// discarding the unrelated Newton, integration, and branch-search knobs.
    pub fn from_period_recovery_config(config: PeriodRecoveryConfig) -> Self {
        Self::new(config.tolerance(), config.agm_max_iterations())
            .expect("validated period-recovery config must induce a valid AGM config")
    }

    /// Returns the tolerance policy used for AGM convergence checks.
    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    /// Returns the AGM iteration budget.
    pub fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    /// Returns the baseline preset for educational AGM experiments.
    pub fn educational_default() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::educational_default())
    }

    /// Returns a tighter preset for more delicate AGM experiments.
    pub fn strict() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::strict())
    }

    /// Returns a more permissive preset for coarse exploratory AGM work.
    pub fn loose() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::loose())
    }
}

/// Terminal outcome for one raw complex AGM run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComplexAgmStatus {
    /// The iterates became close enough under the configured tolerance.
    Succeeded,
    /// The iteration budget was exhausted before the iterates became close.
    HitIterationLimit,
}

/// Which sign of the square root was selected at one AGM step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComplexAgmBranchChoice {
    /// The principal square root of `a_n b_n` was used.
    PrincipalSqrt,
    /// The negative of the principal square root of `a_n b_n` was used.
    NegatedPrincipalSqrt,
}

/// One recorded step of the raw complex AGM iteration.
///
/// Starting from `(a_n, b_n)`, one step forms
///
/// - `next_a = (a_n + b_n) / 2`
/// - `principal_sqrt_product = sqrt(a_n b_n)` on the principal branch
/// - `next_b = ± principal_sqrt_product`
///
/// The sign is chosen deterministically to minimize `|next_a - next_b|`.
/// In exact ties, the principal branch is preferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexAgmIteration {
    index: usize,
    a_n: Complex64,
    b_n: Complex64,
    principal_sqrt_product: Complex64,
    selected_branch: ComplexAgmBranchChoice,
    selected_geometric_mean: Complex64,
    next_a: Complex64,
    next_b: Complex64,
    next_gap_norm: f64,
}

impl ComplexAgmIteration {
    /// Returns the zero-based iteration index.
    pub fn index(&self) -> usize {
        self.index
    }

    /// Returns the arithmetic input `a_n` of this step.
    pub fn a_n(&self) -> &Complex64 {
        &self.a_n
    }

    /// Returns the geometric input `b_n` of this step.
    pub fn b_n(&self) -> &Complex64 {
        &self.b_n
    }

    /// Returns the principal square root `sqrt(a_n b_n)` used as the branch
    /// reference for this step.
    pub fn principal_sqrt_product(&self) -> &Complex64 {
        &self.principal_sqrt_product
    }

    /// Returns which sign branch was selected for the geometric mean.
    pub fn selected_branch(&self) -> ComplexAgmBranchChoice {
        self.selected_branch
    }

    /// Returns the selected signed geometric mean `next_b`.
    pub fn selected_geometric_mean(&self) -> &Complex64 {
        &self.selected_geometric_mean
    }

    /// Returns the next arithmetic mean `(a_n + b_n) / 2`.
    pub fn next_a(&self) -> &Complex64 {
        &self.next_a
    }

    /// Returns the next geometric mean after the branch choice.
    pub fn next_b(&self) -> &Complex64 {
        &self.next_b
    }

    /// Returns the norm `|next_a - next_b|` after this step.
    pub fn next_gap_norm(&self) -> f64 {
        self.next_gap_norm
    }
}

/// Final summary of one raw complex AGM run.
///
/// If the inputs are already close under the configured tolerance, the run can
/// succeed with zero recorded iterations. The reported AGM value is always the
/// symmetric midpoint `(final_a + final_b) / 2` of the last iterate pair.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexAgmResult {
    input_a: Complex64,
    input_b: Complex64,
    agm: Complex64,
    status: ComplexAgmStatus,
    iterations_used: usize,
    final_a: Complex64,
    final_b: Complex64,
    final_gap_norm: f64,
}

impl ComplexAgmResult {
    /// Returns the initial `a` input supplied by the caller.
    pub fn input_a(&self) -> &Complex64 {
        &self.input_a
    }

    /// Returns the initial `b` input supplied by the caller.
    pub fn input_b(&self) -> &Complex64 {
        &self.input_b
    }

    /// Returns the final symmetric AGM approximation.
    pub fn agm(&self) -> &Complex64 {
        &self.agm
    }

    /// Returns how the run terminated.
    pub fn status(&self) -> ComplexAgmStatus {
        self.status
    }

    /// Returns the number of AGM steps actually executed.
    pub fn iterations_used(&self) -> usize {
        self.iterations_used
    }

    /// Returns the final `a_n`.
    pub fn final_a(&self) -> &Complex64 {
        &self.final_a
    }

    /// Returns the final `b_n`.
    pub fn final_b(&self) -> &Complex64 {
        &self.final_b
    }

    /// Returns the final gap norm `|final_a - final_b|`.
    pub fn final_gap_norm(&self) -> f64 {
        self.final_gap_norm
    }

    /// Returns whether the run ended in the succeeded state.
    pub fn succeeded(&self) -> bool {
        self.status == ComplexAgmStatus::Succeeded
    }
}

/// Full educational trace of one raw complex AGM run.
///
/// This stores the validated config, every recorded branch choice, and the
/// final result bundle. Higher-level layers can use this to explain why a
/// particular complex square-root branch was preferred at each step.
#[derive(Clone, Debug, PartialEq)]
pub struct ComplexAgmTrace {
    config: ComplexAgmConfig,
    iterations: Vec<ComplexAgmIteration>,
    result: ComplexAgmResult,
}

impl ComplexAgmTrace {
    /// Returns the validated config used for this trace.
    pub fn config(&self) -> ComplexAgmConfig {
        self.config
    }

    /// Returns the recorded AGM steps in execution order.
    pub fn iterations(&self) -> &[ComplexAgmIteration] {
        &self.iterations
    }

    /// Returns the final result bundle.
    pub fn result(&self) -> &ComplexAgmResult {
        &self.result
    }
}

/// Runs the raw complex AGM iteration without storing the full step trace.
///
/// For inputs `(a, b)`, this routine repeatedly applies
/// `a_{n+1} = (a_n + b_n) / 2` and chooses `b_{n+1} = ± sqrt(a_n b_n)`
/// so that `|a_{n+1} - b_{n+1}|` is minimized. The principal square root is
/// used as the reference branch, and exact ties prefer that principal branch.
///
/// The run succeeds once `a_n` and `b_n` are close under the mixed
/// absolute/relative comparison policy stored in `config`.
///
/// Complexity: `Θ(n)`, where `n` is the maximum number of iterations.
pub fn complex_agm(
    a: Complex64,
    b: Complex64,
    config: ComplexAgmConfig,
) -> Result<ComplexAgmResult, AnalyticCurveError> {
    run_complex_agm(a, b, config, false).map(|(_, result)| result)
}

/// Runs the raw complex AGM iteration and records every step.
///
/// This is the pedagogical companion to [`complex_agm`]. It stores the
/// principal square root, selected sign branch, and post-step gap norm at each
/// iteration so later higher-level routines can explain branch choices.
///
/// Complexity: `Θ(n)`, where `n` is the maximum number of iterations.
pub fn complex_agm_trace(
    a: Complex64,
    b: Complex64,
    config: ComplexAgmConfig,
) -> Result<ComplexAgmTrace, AnalyticCurveError> {
    let (iterations, result) = run_complex_agm(a, b, config, true)?;

    Ok(ComplexAgmTrace {
        config,
        iterations,
        result,
    })
}

fn run_complex_agm(
    a: Complex64,
    b: Complex64,
    config: ComplexAgmConfig,
    record_trace: bool,
) -> Result<(Vec<ComplexAgmIteration>, ComplexAgmResult), AnalyticCurveError> {
    if !a.is_finite() || !b.is_finite() {
        return Err(AnalyticCurveError::InvalidAgmInput);
    }

    let mut current_a = a;
    let mut current_b = b;
    let mut iterations = Vec::new();
    let initial_gap = (current_a - current_b).norm();

    if ComplexApprox::eq_with_tolerance(&current_a, &current_b, config.tolerance()) {
        return Ok((
            iterations,
            ComplexAgmResult {
                input_a: a,
                input_b: b,
                agm: (current_a + current_b) / 2.0,
                status: ComplexAgmStatus::Succeeded,
                iterations_used: 0,
                final_a: current_a,
                final_b: current_b,
                final_gap_norm: initial_gap,
            },
        ));
    }

    for index in 0..config.max_iterations() {
        let step = next_agm_iteration(index, current_a, current_b);
        current_a = step.next_a;
        current_b = step.next_b;

        if record_trace {
            iterations.push(step);
        }

        if ComplexApprox::eq_with_tolerance(&current_a, &current_b, config.tolerance()) {
            return Ok((
                iterations,
                ComplexAgmResult {
                    input_a: a,
                    input_b: b,
                    agm: (current_a + current_b) / 2.0,
                    status: ComplexAgmStatus::Succeeded,
                    iterations_used: index + 1,
                    final_a: current_a,
                    final_b: current_b,
                    final_gap_norm: (current_a - current_b).norm(),
                },
            ));
        }
    }

    Ok((
        iterations,
        ComplexAgmResult {
            input_a: a,
            input_b: b,
            agm: (current_a + current_b) / 2.0,
            status: ComplexAgmStatus::HitIterationLimit,
            iterations_used: config.max_iterations(),
            final_a: current_a,
            final_b: current_b,
            final_gap_norm: (current_a - current_b).norm(),
        },
    ))
}

fn next_agm_iteration(index: usize, a_n: Complex64, b_n: Complex64) -> ComplexAgmIteration {
    let next_a = (a_n + b_n) / 2.0;
    let principal_sqrt_product = (a_n * b_n).sqrt();
    let negated_principal = -principal_sqrt_product;
    let principal_gap = (next_a - principal_sqrt_product).norm();
    let negated_gap = (next_a - negated_principal).norm();

    let (selected_branch, selected_geometric_mean, next_gap_norm) = if principal_gap <= negated_gap
    {
        (
            ComplexAgmBranchChoice::PrincipalSqrt,
            principal_sqrt_product,
            principal_gap,
        )
    } else {
        (
            ComplexAgmBranchChoice::NegatedPrincipalSqrt,
            negated_principal,
            negated_gap,
        )
    };

    ComplexAgmIteration {
        index,
        a_n,
        b_n,
        principal_sqrt_product,
        selected_branch,
        selected_geometric_mean,
        next_a,
        next_b: selected_geometric_mean,
        next_gap_norm,
    }
}
