use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::{CharacterSumPointCount, FrobeniusTrace, HasseInterval};
use crate::elliptic_curves::short_weierstrass::{
    ExponentLowerBoundGroupOrderVerification, PointOrderFromMultipleReport,
};
use crate::fields::FiniteFieldDescriptor;
use num_bigint::BigUint;

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
pub enum MestreSide {
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
pub struct MestreStepReport {
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
    pub fn side(&self) -> MestreSide {
        self.side
    }

    /// Returns the annihilating multiple found in `H(p)` for the sampled
    /// point.
    pub fn annihilating_multiple(&self) -> u128 {
        self.annihilating_multiple
    }

    /// Returns the exact point-order report recovered from that multiple.
    pub fn point_order_report(&self) -> &PointOrderFromMultipleReport {
        &self.point_order_report
    }

    /// Returns the running exponent lower bound after this step.
    pub fn accumulated_exponent_lower_bound(&self) -> &BigUint {
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

    /// Returns which side first produced the decisive unique multiple in
    /// `H(p)`.
    pub fn resolved_side(&self) -> MestreSide {
        self.resolved_side
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

    /// Returns the recorded alternating Mestre steps.
    pub fn steps(&self) -> &[MestreStepReport] {
        &self.steps
    }
}

/// Public strategy choices for counting `#E(F_q)` on small finite curves.
///
/// The current educational implementation distinguishes:
///
/// - [`Self::Exhaustive`], which materializes all rational points directly
/// - [`Self::QuadraticCharacter`], which uses `#E(F_q) = q + 1 + Σ_{x ∈ F_q} χ(f(x))`
/// - [`Self::Auto`], which chooses the best implemented route for the curve
///
/// At the moment [`Self::Auto`] prefers the quadratic-character route when the
/// curve implements it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HasseGroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
}

impl HasseGroupOrderStrategy {
    pub(crate) fn as_group_order_strategy(self) -> GroupOrderStrategy {
        match self {
            Self::Auto => GroupOrderStrategy::Auto,
            Self::Exhaustive => GroupOrderStrategy::Exhaustive,
            Self::QuadraticCharacter => GroupOrderStrategy::QuadraticCharacter,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupOrderStrategy {
    Auto,
    Exhaustive,
    QuadraticCharacter,
    /// Prime-field Mestre route from Lecture 7, Algorithm 7.8.
    ///
    /// This strategy is specific to curves over `F_p` and alternates between
    /// the original curve and one quadratic twist while accumulating lower
    /// bounds for their exponents until one side has a unique multiple in the
    /// Hasse interval.
    MestreFp(MestreConfig),
    /// Tries to recover `#E(F_q)` from one known lower bound for the exponent
    /// `λ(E(F_q))` together with a Hasse-interval-producing base strategy.
    ///
    /// This route does not recount points from scratch. Instead it first runs
    /// the chosen `hasse_strategy` to obtain one Hasse interval `H(q)`, then
    /// asks whether that interval contains exactly one multiple of the
    /// supplied lower bound. It succeeds only in that uniqueness case.
    FromExponentLowerBoundAndPointCount {
        exponent_lower_bound: BigUint,
        hasse_strategy: HasseGroupOrderStrategy,
    },
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
    MestreFp(Box<MestreGroupOrderReport>),
    FromExponentLowerBound(Box<ExponentLowerBoundGroupOrderVerification>),
}

impl GroupOrderReport {
    /// Returns the strategy used to build this report.
    pub fn strategy(&self) -> GroupOrderStrategy {
        match self {
            Self::ExhaustiveTrace(_) => GroupOrderStrategy::Exhaustive,
            Self::QuadraticCharacter(_) => GroupOrderStrategy::QuadraticCharacter,
            Self::MestreFp(report) => GroupOrderStrategy::MestreFp(report.config().clone()),
            Self::FromExponentLowerBound(report) => {
                GroupOrderStrategy::FromExponentLowerBoundAndPointCount {
                    exponent_lower_bound: report.exponent_lower_bound().clone(),
                    hasse_strategy: match report.group_order_report().strategy() {
                        GroupOrderStrategy::Auto => HasseGroupOrderStrategy::Auto,
                        GroupOrderStrategy::Exhaustive => HasseGroupOrderStrategy::Exhaustive,
                        GroupOrderStrategy::QuadraticCharacter => {
                            HasseGroupOrderStrategy::QuadraticCharacter
                        }
                        GroupOrderStrategy::MestreFp(_) => {
                            unreachable!(
                                "the lower-bound report should wrap one Hasse-based strategy"
                            )
                        }
                        GroupOrderStrategy::FromExponentLowerBoundAndPointCount { .. } => {
                            unreachable!(
                                "the lower-bound report should wrap one base group-order strategy"
                            )
                        }
                    },
                }
            }
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => trace.field_order(),
            Self::QuadraticCharacter(report) => report.field_order(),
            Self::MestreFp(report) => report.field_order(),
            Self::FromExponentLowerBound(report) => report.group_order_report().field_order(),
        }
    }

    /// Returns the counted curve order `#E(F_q)`.
    pub fn curve_order(&self) -> u128 {
        match self {
            Self::ExhaustiveTrace(trace) => u128::from(trace.curve_order()),
            Self::QuadraticCharacter(report) => report.curve_order(),
            Self::MestreFp(report) => report.curve_order(),
            Self::FromExponentLowerBound(report) => report.group_order_report().curve_order(),
        }
    }

    /// Returns the Frobenius trace `t = q + 1 - #E(F_q)`.
    pub fn trace(&self) -> i128 {
        match self {
            Self::ExhaustiveTrace(trace) => i128::from(trace.trace()),
            Self::QuadraticCharacter(report) => report.trace(),
            Self::MestreFp(report) => report.trace(),
            Self::FromExponentLowerBound(report) => report.group_order_report().trace(),
        }
    }

    /// Returns the discrete Hasse interval attached to the same `F_q`.
    pub fn hasse_interval(&self) -> HasseInterval {
        match self {
            Self::ExhaustiveTrace(trace) => trace.hasse_interval(),
            Self::QuadraticCharacter(report) => report.hasse_interval(),
            Self::MestreFp(report) => report.hasse_interval(),
            Self::FromExponentLowerBound(report) => report.group_order_report().hasse_interval(),
        }
    }

    /// Converts this report into the shared Frobenius-trace package.
    pub fn to_frobenius_trace(&self) -> Result<FrobeniusTrace, CurveError> {
        match self {
            Self::ExhaustiveTrace(trace) => Ok(trace.clone()),
            Self::QuadraticCharacter(report) => report.to_frobenius_trace(),
            Self::MestreFp(report) => Ok(report.original().clone()),
            Self::FromExponentLowerBound(report) => {
                report.group_order_report().to_frobenius_trace()
            }
        }
    }
}
