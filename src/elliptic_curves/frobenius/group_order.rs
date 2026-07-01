use num_bigint::BigUint;

use crate::elliptic_curves::{
    CurveError,
    frobenius::schoof::{SchoofGroupOrderOutcome, SchoofGroupOrderReport},
    frobenius::{FrobeniusTrace, HasseInterval, character_sum::CharacterSumPointCount},
    short_weierstrass::point_order::PointOrderFromMultipleReport,
};
use crate::fields::{finite_field_descriptor::FiniteFieldDescriptor, traits::FiniteField};

/// Configuration for the prime-field Mestre group-order route.
///
/// It is a Las Vegas algorithm over `F_p` that alternates between a curve and
/// one quadratic twist while accumulating lower bounds for the corresponding
/// exponents.
///
/// The current surface exposes only one optional iteration cap:
/// - `None` means "keep sampling until the route certifies one unique order"
/// - `Some(k)` means "abort after at most `k` Mestre iterations"
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MestreConfig {
    max_iterations: Option<usize>,
}

impl MestreConfig {
    /// Returns a config with no iteration cap.
    pub const fn unbounded() -> Self {
        Self {
            max_iterations: None,
        }
    }

    /// Returns a config that stops after at most `max_iterations` steps.
    pub const fn with_iteration_cap(max_iterations: usize) -> Self {
        Self {
            max_iterations: Some(max_iterations),
        }
    }

    /// Returns the optional cap on Mestre iterations.
    pub const fn max_iterations(&self) -> Option<usize> {
        self.max_iterations
    }
}

/// Which side of Mestre's theorem produced the decisive unique multiple in
/// the Hasse interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MestreSide {
    Original,
    QuadraticTwist,
}

/// One alternating-step record in Mestre's algorithm.
///
/// Each step samples a point on either the original curve or its quadratic
/// twist, finds an annihilating multiple in `H(p)`, recovers the exact point
/// order from that multiple, and updates the running exponent lower bound on
/// that side.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MestreStepReport {
    side: MestreSide,
    annihilating_multiple: u128,
    point_order_report: PointOrderFromMultipleReport,
    accumulated_exponent_lower_bound: BigUint,
}

impl MestreStepReport {
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn new(
        side: MestreSide,
        annihilating_multiple: u128,
        point_order_report: PointOrderFromMultipleReport,
        accumulated_exponent_lower_bound: BigUint,
    ) -> Self {
        Self {
            side,
            annihilating_multiple,
            point_order_report,
            accumulated_exponent_lower_bound,
        }
    }

    /// Returns whether this step ran on the original curve or on its twist.
    pub(crate) fn side(&self) -> MestreSide {
        self.side
    }

    /// Returns the annihilating multiple found in `H(p)` for the sampled
    /// point.
    #[cfg(feature = "visualization")]
    pub(crate) fn annihilating_multiple(&self) -> u128 {
        self.annihilating_multiple
    }

    /// Returns the exact point-order report recovered from that multiple.
    #[cfg(feature = "visualization")]
    pub(crate) fn point_order_report(&self) -> &PointOrderFromMultipleReport {
        &self.point_order_report
    }

    /// Returns the running exponent lower bound after this step.
    pub(crate) fn accumulated_exponent_lower_bound(&self) -> &BigUint {
        &self.accumulated_exponent_lower_bound
    }
}

/// Route-preserving report for the prime-field Mestre group-order algorithm.
///
/// This report keeps both sides of the algorithm visible:
/// - the Frobenius package for the original curve `E/F_p`
/// - the Frobenius package for the chosen quadratic twist `E'/F_p`
/// - the alternating step history
///
/// The final resolved group order returned through the surrounding
/// [`GroupOrderReport`] is always `#E(F_p)` for the original curve, even when
/// uniqueness was first certified on the twist side. The running exponent
/// lower bounds are derived from the last recorded step on each side, so
/// `steps` is the single source of truth for Mestre's iterative progress.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MestreGroupOrderReport {
    config: MestreConfig,
    resolved_side: MestreSide,
    original: FrobeniusTrace,
    twist: FrobeniusTrace,
    steps: Vec<MestreStepReport>,
}

impl MestreGroupOrderReport {
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn new(
        config: MestreConfig,
        resolved_side: MestreSide,
        original: FrobeniusTrace,
        twist: FrobeniusTrace,
        steps: Vec<MestreStepReport>,
    ) -> Self {
        debug_assert_eq!(original.base_field(), twist.base_field());

        Self {
            config,
            resolved_side,
            original,
            twist,
            steps,
        }
    }

    fn exponent_lower_bound_for_side(&self, side: MestreSide) -> BigUint {
        self.steps
            .iter()
            .rev()
            .find(|step| step.side() == side)
            .map(|step| step.accumulated_exponent_lower_bound().clone())
            .unwrap_or_else(|| BigUint::from(1u8))
    }

    /// Returns the configuration used by the Mestre route.
    pub fn config(&self) -> &MestreConfig {
        &self.config
    }

    /// Returns the Frobenius-trace package for the original curve `E/F_p`.
    pub fn original(&self) -> &FrobeniusTrace {
        &self.original
    }

    /// Returns the Frobenius-trace package for the chosen quadratic twist.
    pub fn twist(&self) -> &FrobeniusTrace {
        &self.twist
    }

    /// Returns the common finite base-field descriptor.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        self.original.base_field()
    }

    /// Returns the base-field order `p`.
    pub fn field_order(&self) -> u128 {
        self.original.field_order()
    }

    /// Returns `#E(F_p)` for the original curve.
    pub fn curve_order(&self) -> u128 {
        u128::from(self.original.curve_order())
    }

    /// Returns `#E'(F_p)` for the chosen quadratic twist.
    pub fn twist_curve_order(&self) -> u128 {
        u128::from(self.twist.curve_order())
    }

    /// Returns the Frobenius trace `t = p + 1 - #E(F_p)` of the original
    /// curve.
    pub fn trace(&self) -> i128 {
        i128::from(self.original.trace())
    }

    /// Returns the shared Hasse interval `H(p)`.
    pub fn hasse_interval(&self) -> HasseInterval {
        self.original.hasse_interval()
    }

    /// Returns the accumulated lower bound for `λ(E(F_p))`.
    pub fn original_exponent_lower_bound(&self) -> BigUint {
        self.exponent_lower_bound_for_side(MestreSide::Original)
    }

    /// Returns the accumulated lower bound for `λ(E'(F_p))`.
    pub fn twist_exponent_lower_bound(&self) -> BigUint {
        self.exponent_lower_bound_for_side(MestreSide::QuadraticTwist)
    }

    /// Returns which side first produced the decisive unique multiple in
    /// `H(p)`, as a stable human-facing label.
    pub fn resolved_side_label(&self) -> &'static str {
        match self.resolved_side {
            MestreSide::Original => "original curve",
            MestreSide::QuadraticTwist => "quadratic twist",
        }
    }

    /// Returns the group-order candidate on the side that first became unique
    /// in `H(p)`.
    ///
    /// This equals `#E(F_p)` when the original curve resolved first and
    /// `#E'(F_p)` when the quadratic twist resolved first.
    pub fn resolved_side_group_order_candidate(&self) -> u128 {
        match self.resolved_side {
            MestreSide::Original => self.curve_order(),
            MestreSide::QuadraticTwist => self.twist_curve_order(),
        }
    }

    /// Returns how many alternating Mestre steps were recorded.
    pub fn step_count(&self) -> usize {
        self.steps.len()
    }

    /// Returns the recorded alternating Mestre steps.
    #[cfg(feature = "visualization")]
    pub(crate) fn steps(&self) -> &[MestreStepReport] {
        &self.steps
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteFieldGroupOrderStrategy {
    Auto,
    Schoof,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmallFieldGroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
    Schoof,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SmallFieldSampledGroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
    Schoof,
    /// This strategy is specific to curves over `F_p` and alternates between
    /// the original curve and one quadratic twist while accumulating lower
    /// bounds for their exponents until one side has a unique multiple in the
    /// Hasse interval.
    MestreFp(MestreConfig),
}

/// Route label attached to one realized [`GroupOrderReport`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GroupOrderRoute {
    Exhaustive,
    QuadraticCharacter,
    Schoof,
    MestreFp,
}

/// Route-preserving summary for the automatic Schoof group-order strategy.
///
/// The fully detailed Schoof arithmetic report stays available under
/// `elliptic_curves::frobenius::schoof`, where it remains generic in the base
/// field and preserves quotient-ring data. This summary is the non-generic
/// surface used by the shared `GroupOrderReport` enum.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SchoofGroupOrderSummary {
    resolved: FrobeniusTrace,
    attempted_odd_primes: Vec<usize>,
    combined_crt_modulus: BigUint,
}

impl SchoofGroupOrderSummary {
    pub(crate) fn from_detailed<F: FiniteField>(
        report: &SchoofGroupOrderReport<F>,
    ) -> Result<Self, CurveError> {
        match report.outcome() {
            SchoofGroupOrderOutcome::GroupOrderFound { .. } => {
                let resolved = report
                    .to_frobenius_trace()?
                    .expect("successful Schoof group-order reports should carry a trace");
                let attempted_odd_primes = report
                    .crt_report()
                    .odd_prime_reports()
                    .iter()
                    .map(|odd_prime_report| odd_prime_report.odd_prime())
                    .collect();
                let combined_crt_modulus = report
                    .crt_report()
                    .combined_solution()
                    .expect("successful Schoof group-order reports should carry one CRT solution")
                    .modulus()
                    .clone();
                Ok(Self {
                    resolved,
                    attempted_odd_primes,
                    combined_crt_modulus,
                })
            }
            SchoofGroupOrderOutcome::BlockedOnOddPrime => {
                let blocked_prime = report
                    .crt_report()
                    .odd_prime_reports()
                    .last()
                    .map(|odd_prime_report| odd_prime_report.odd_prime())
                    .expect(
                        "blocked Schoof route should have recorded the blocking odd-prime step",
                    );
                Err(CurveError::SchoofBlockedOnOddPrime {
                    odd_prime: blocked_prime,
                })
            }
            SchoofGroupOrderOutcome::AmbiguousTraceClass {
                candidate_count, ..
            } => Err(CurveError::SchoofAmbiguousTraceClass {
                candidate_count: *candidate_count,
            }),
            SchoofGroupOrderOutcome::InconsistentWithHasse => {
                Err(CurveError::SchoofInconsistentWithHasse)
            }
        }
    }

    /// Returns the resolved Frobenius-trace package.
    pub fn resolved(&self) -> &FrobeniusTrace {
        &self.resolved
    }

    /// Returns the odd primes attempted before Schoof resolved the trace.
    pub fn attempted_odd_primes(&self) -> &[usize] {
        &self.attempted_odd_primes
    }

    /// Returns the final CRT modulus that already exceeded the Hasse
    /// uniqueness threshold.
    pub fn combined_crt_modulus(&self) -> &BigUint {
        &self.combined_crt_modulus
    }
}

/// Shared group-order result returned by curve-side methods.
///
/// This enum keeps the algorithmic route explicit instead of collapsing every
/// count immediately into a bare integer. Callers can still recover the common
/// arithmetic data such as `#E(F_q)`, the trace `t`, and the Hasse interval
/// without losing which strategy produced it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupOrderReport {
    ExhaustiveTrace(FrobeniusTrace),
    QuadraticCharacter(CharacterSumPointCount),
    Schoof(Box<SchoofGroupOrderSummary>),
    MestreFp(Box<MestreGroupOrderReport>),
}

impl GroupOrderReport {
    /// Returns the algorithmic route that produced this report.
    pub fn route(&self) -> GroupOrderRoute {
        match self {
            Self::ExhaustiveTrace(_) => GroupOrderRoute::Exhaustive,
            Self::QuadraticCharacter(_) => GroupOrderRoute::QuadraticCharacter,
            Self::Schoof(_) => GroupOrderRoute::Schoof,
            Self::MestreFp(_) => GroupOrderRoute::MestreFp,
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => trace.field_order(),
            Self::QuadraticCharacter(report) => report.field_order(),
            Self::Schoof(report) => report.resolved().field_order(),
            Self::MestreFp(report) => report.field_order(),
        }
    }

    /// Returns the counted curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => u128::from(trace.curve_order()),
            Self::QuadraticCharacter(report) => report.curve_order(),
            Self::Schoof(report) => u128::from(report.resolved().curve_order()),
            Self::MestreFp(report) => report.curve_order(),
        }
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> i128 {
        match self {
            Self::ExhaustiveTrace(trace) => i128::from(trace.trace()),
            Self::QuadraticCharacter(report) => report.trace(),
            Self::Schoof(report) => i128::from(report.resolved().trace()),
            Self::MestreFp(report) => report.trace(),
        }
    }

    /// Returns the discrete Hasse interval attached to the same `F_q`.
    pub fn hasse_interval(&self) -> HasseInterval {
        match self {
            Self::ExhaustiveTrace(trace) => trace.hasse_interval(),
            Self::QuadraticCharacter(report) => report.hasse_interval(),
            Self::Schoof(report) => report.resolved().hasse_interval(),
            Self::MestreFp(report) => report.hasse_interval(),
        }
    }

    /// Converts this report into the shared Frobenius-trace package.
    pub fn to_frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        match self {
            Self::ExhaustiveTrace(trace) => Ok(trace.clone()),
            Self::QuadraticCharacter(report) => report.to_frobenius_trace(),
            Self::Schoof(report) => Ok(report.resolved().clone()),
            Self::MestreFp(report) => Ok(report.original().clone()),
        }
    }
}
