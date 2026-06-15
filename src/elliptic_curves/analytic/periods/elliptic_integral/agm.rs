use num_complex::Complex64;

use crate::elliptic_curves::analytic::{AnalyticCurveError, periods::config::PeriodRecoveryConfig};
use crate::fields::complex_approx::ComplexApprox;
use crate::numerics::ApproxTolerance;

/// Validated numerical policy for the raw complex AGM iteration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexAgmConfig {
    tolerance: ApproxTolerance,
    max_iterations: usize,
}

impl ComplexAgmConfig {
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

    pub fn from_period_recovery_config(config: PeriodRecoveryConfig) -> Self {
        Self::new(config.tolerance(), config.agm_max_iterations())
            .expect("validated period-recovery config must induce a valid AGM config")
    }

    pub fn tolerance(&self) -> ApproxTolerance {
        self.tolerance
    }

    pub fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    pub fn educational_default() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::educational_default())
    }

    pub fn strict() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::strict())
    }

    pub fn loose() -> Self {
        Self::from_period_recovery_config(PeriodRecoveryConfig::loose())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComplexAgmStatus {
    Succeeded,
    HitIterationLimit,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComplexAgmBranchChoice {
    PrincipalSqrt,
    NegatedPrincipalSqrt,
}

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
    fn from_state(index: usize, a_n: Complex64, b_n: Complex64) -> Self {
        let next_a = (a_n + b_n) / 2.0;
        let principal_sqrt_product = (a_n * b_n).sqrt();
        let negated_principal = -principal_sqrt_product;
        let principal_gap = (next_a - principal_sqrt_product).norm();
        let negated_gap = (next_a - negated_principal).norm();

        let (selected_branch, selected_geometric_mean, next_gap_norm) =
            if principal_gap <= negated_gap {
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

        Self {
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

    pub fn index(&self) -> usize {
        self.index
    }
    pub fn a_n(&self) -> &Complex64 {
        &self.a_n
    }
    pub fn b_n(&self) -> &Complex64 {
        &self.b_n
    }
    pub fn principal_sqrt_product(&self) -> &Complex64 {
        &self.principal_sqrt_product
    }
    pub fn selected_branch(&self) -> ComplexAgmBranchChoice {
        self.selected_branch
    }
    pub fn selected_geometric_mean(&self) -> &Complex64 {
        &self.selected_geometric_mean
    }
    pub fn next_a(&self) -> &Complex64 {
        &self.next_a
    }
    pub fn next_b(&self) -> &Complex64 {
        &self.next_b
    }
    pub fn next_gap_norm(&self) -> f64 {
        self.next_gap_norm
    }
}

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
    pub fn input_a(&self) -> &Complex64 {
        &self.input_a
    }
    pub fn input_b(&self) -> &Complex64 {
        &self.input_b
    }
    pub fn agm(&self) -> &Complex64 {
        &self.agm
    }
    pub fn status(&self) -> ComplexAgmStatus {
        self.status
    }
    pub fn iterations_used(&self) -> usize {
        self.iterations_used
    }
    pub fn final_a(&self) -> &Complex64 {
        &self.final_a
    }
    pub fn final_b(&self) -> &Complex64 {
        &self.final_b
    }
    pub fn final_gap_norm(&self) -> f64 {
        self.final_gap_norm
    }
    pub fn succeeded(&self) -> bool {
        self.status == ComplexAgmStatus::Succeeded
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ComplexAgmTrace {
    config: ComplexAgmConfig,
    iterations: Vec<ComplexAgmIteration>,
    result: ComplexAgmResult,
}

impl ComplexAgmTrace {
    fn run(
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
            let step = ComplexAgmIteration::from_state(index, current_a, current_b);
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

    pub fn config(&self) -> ComplexAgmConfig {
        self.config
    }
    pub fn iterations(&self) -> &[ComplexAgmIteration] {
        &self.iterations
    }
    pub fn result(&self) -> &ComplexAgmResult {
        &self.result
    }
}

pub fn complex_agm(
    a: Complex64,
    b: Complex64,
    config: ComplexAgmConfig,
) -> Result<ComplexAgmResult, AnalyticCurveError> {
    ComplexAgmTrace::run(a, b, config, false).map(|(_, result)| result)
}

pub fn complex_agm_trace(
    a: Complex64,
    b: Complex64,
    config: ComplexAgmConfig,
) -> Result<ComplexAgmTrace, AnalyticCurveError> {
    let (iterations, result) = ComplexAgmTrace::run(a, b, config, true)?;
    Ok(ComplexAgmTrace {
        config,
        iterations,
        result,
    })
}
