use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::{
    CurveError,
    frobenius::schoof::{SchoofTraceCrtOutcome, SchoofTraceCrtReport},
    frobenius::{FrobeniusTrace, HasseInterval},
};
use crate::fields::finite_field_descriptor::FiniteFieldDescriptor;
use crate::fields::traits::FiniteField;
use crate::numerics::{ceil_div_bigint_by_positive, floor_div_bigint_by_positive};

/// Final outcome of the current Schoof group-order stage.
///
/// This stage starts from one CRT class for the Frobenius trace `t` and then
/// intersects that class with Hasse's trace interval `[-⌊2√q⌋, ⌊2√q⌋]`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SchoofGroupOrderOutcome {
    /// The CRT class meets Hasse's trace interval in exactly one integer.
    GroupOrderFound { trace: BigInt, curve_order: BigUint },
    /// The CRT class meets Hasse's trace interval in several integers, so
    /// more odd primes are still needed.
    AmbiguousTraceClass {
        first_trace: BigInt,
        last_trace: BigInt,
        candidate_count: BigUint,
    },
    /// The currently available CRT stage stopped before producing one full
    /// trace class modulo the product of the requested primes.
    BlockedOnOddPrime,
    /// The CRT class did not meet Hasse's trace interval at all.
    ///
    /// Mathematically this should not happen for a correct Schoof run.
    InconsistentWithHasse,
}

/// Route-preserving report for the current end-to-end Schoof group-order stage.
///
/// The stored CRT report remains the source of truth for the prime-by-prime
/// Schoof steps. This wrapper adds only the final Hasse-resolution summary
/// needed to recover `t` and then `#E(F_q) = q + 1 - t`.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofGroupOrderReport<F: FiniteField> {
    base_field: FiniteFieldDescriptor,
    crt_report: SchoofTraceCrtReport<F>,
    outcome: SchoofGroupOrderOutcome,
}

impl<F: FiniteField> SchoofGroupOrderReport<F> {
    pub(crate) fn new(
        base_field: FiniteFieldDescriptor,
        crt_report: SchoofTraceCrtReport<F>,
        outcome: SchoofGroupOrderOutcome,
    ) -> Self {
        Self {
            base_field,
            crt_report,
            outcome,
        }
    }

    /// Returns the finite base-field descriptor.
    pub fn base_field(&self) -> &FiniteFieldDescriptor {
        &self.base_field
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> BigUint {
        self.crt_report.field_order().clone()
    }

    /// Returns the underlying CRT-stage Schoof report.
    pub fn crt_report(&self) -> &SchoofTraceCrtReport<F> {
        &self.crt_report
    }

    /// Returns the discrete Hasse interval `H(q)` for `#E(F_q)`.
    ///
    /// Complexity: `Θ(1)`.
    pub fn hasse_interval(&self) -> HasseInterval {
        HasseInterval::for_q(self.crt_report.field_order().clone())
            .expect("stored Schoof field order should define a valid Hasse interval")
    }

    /// Returns the final Hasse-resolution outcome.
    pub fn outcome(&self) -> &SchoofGroupOrderOutcome {
        &self.outcome
    }

    /// Returns the recovered curve order `#E(F_q)` when the current Schoof
    /// data already determines it uniquely.
    pub fn curve_order(&self) -> Option<&BigUint> {
        match &self.outcome {
            SchoofGroupOrderOutcome::GroupOrderFound { curve_order, .. } => Some(curve_order),
            _ => None,
        }
    }

    /// Returns the recovered Frobenius trace `t = q + 1 - #E(F_q)` when the
    /// current Schoof data already determines it uniquely.
    pub fn trace(&self) -> Option<&BigInt> {
        match &self.outcome {
            SchoofGroupOrderOutcome::GroupOrderFound { trace, .. } => Some(trace),
            _ => None,
        }
    }

    /// Converts the resolved Schoof output into the shared `FrobeniusTrace`
    /// package when the current data already determines one unique trace.
    ///
    /// Complexity: `Θ(1)`.
    pub fn to_frobenius_trace(&self) -> Result<Option<FrobeniusTrace>, CurveError> {
        let Some(curve_order) = self.curve_order() else {
            return Ok(None);
        };
        FrobeniusTrace::from_order(self.base_field.clone(), curve_order.clone()).map(Some)
    }
}

pub(crate) fn finalize_schoof_group_order_report<F: FiniteField>(
    base_field: FiniteFieldDescriptor,
    crt_report: SchoofTraceCrtReport<F>,
) -> Result<SchoofGroupOrderReport<F>, CurveError> {
    let outcome = match crt_report.outcome() {
        SchoofTraceCrtOutcome::BlockedOnOddPrime { .. } => {
            SchoofGroupOrderOutcome::BlockedOnOddPrime
        }
        SchoofTraceCrtOutcome::Combined { solution } => {
            let compatible = hasse_compatible_trace_class(crt_report.field_order(), solution)?;
            match compatible {
                HasseCompatibleTraceClass::Unique { trace } => {
                    let field_order = crt_report.field_order().clone();
                    let curve_order = crate::elliptic_curves::frobenius::invariants::curve_order_from_field_order_and_trace(
                            &field_order,
                            &trace,
                        )?;
                    SchoofGroupOrderOutcome::GroupOrderFound { trace, curve_order }
                }
                HasseCompatibleTraceClass::Ambiguous {
                    first_trace,
                    last_trace,
                    candidate_count,
                } => SchoofGroupOrderOutcome::AmbiguousTraceClass {
                    first_trace,
                    last_trace,
                    candidate_count,
                },
                HasseCompatibleTraceClass::Empty => SchoofGroupOrderOutcome::InconsistentWithHasse,
            }
        }
    };

    Ok(SchoofGroupOrderReport::new(base_field, crt_report, outcome))
}

enum HasseCompatibleTraceClass {
    Empty,
    Unique {
        trace: BigInt,
    },
    Ambiguous {
        first_trace: BigInt,
        last_trace: BigInt,
        candidate_count: BigUint,
    },
}

fn hasse_compatible_trace_class(
    field_order: &BigUint,
    solution: &crate::numerics::chinese_remainder::ChineseRemainderSolution,
) -> Result<HasseCompatibleTraceClass, CurveError> {
    let hasse_interval = HasseInterval::for_q(field_order.clone())?;
    let trace_bound = BigInt::from(hasse_interval.trace_bound());
    let lower = -trace_bound.clone();
    let upper = trace_bound;
    let residue = BigInt::from(solution.residue().clone());
    let modulus = BigInt::from(solution.modulus().clone());

    let k_min = ceil_div_bigint_by_positive(&(lower - &residue), &modulus);
    let k_max = floor_div_bigint_by_positive(&(upper - &residue), &modulus);
    if k_min > k_max {
        return Ok(HasseCompatibleTraceClass::Empty);
    }

    let first_trace_bigint = &residue + (&k_min * &modulus);
    let last_trace_bigint = &residue + (&k_max * &modulus);
    let candidate_count = (&k_max - &k_min + BigInt::from(1u8)).to_biguint().ok_or(
        CurveError::InvalidHasseIntervalFieldOrder {
            field_order: field_order.clone(),
        },
    )?;

    if candidate_count == BigUint::from(1u8) {
        Ok(HasseCompatibleTraceClass::Unique {
            trace: first_trace_bigint,
        })
    } else {
        Ok(HasseCompatibleTraceClass::Ambiguous {
            first_trace: first_trace_bigint,
            last_trace: last_trace_bigint,
            candidate_count,
        })
    }
}
