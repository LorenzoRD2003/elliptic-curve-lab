use crate::elliptic_curves::frobenius::FrobeniusTrace;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicEquationTerms<P> {
    pi_q: P,
    pi_q_squared: P,
    trace_term: P,
    q_times_point: P,
    lhs: P,
}

impl<P> FrobeniusCharacteristicEquationTerms<P> {
    pub(crate) fn new(pi_q: P, pi_q_squared: P, trace_term: P, q_times_point: P, lhs: P) -> Self {
        Self {
            pi_q,
            pi_q_squared,
            trace_term,
            q_times_point,
            lhs,
        }
    }

    pub fn pi_q(&self) -> &P {
        &self.pi_q
    }

    pub fn pi_q_squared(&self) -> &P {
        &self.pi_q_squared
    }

    pub fn trace_term(&self) -> &P {
        &self.trace_term
    }

    pub fn q_times_point(&self) -> &P {
        &self.q_times_point
    }

    pub fn lhs(&self) -> &P {
        &self.lhs
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicEquationCheck<P> {
    point: P,
    pi_q: P,
    pi_q_squared: P,
    trace_term: P,
    q_times_point: P,
    lhs: P,
    holds: bool,
}

impl<P> FrobeniusCharacteristicEquationCheck<P> {
    pub(crate) fn from_terms(
        point: P,
        terms: FrobeniusCharacteristicEquationTerms<P>,
        holds: bool,
    ) -> Self {
        Self {
            point,
            pi_q: terms.pi_q,
            pi_q_squared: terms.pi_q_squared,
            trace_term: terms.trace_term,
            q_times_point: terms.q_times_point,
            lhs: terms.lhs,
            holds,
        }
    }

    pub fn point(&self) -> &P {
        &self.point
    }

    pub fn pi_q(&self) -> &P {
        &self.pi_q
    }

    pub fn pi_q_squared(&self) -> &P {
        &self.pi_q_squared
    }

    pub fn trace_term(&self) -> &P {
        &self.trace_term
    }

    pub fn q_times_point(&self) -> &P {
        &self.q_times_point
    }

    pub fn lhs(&self) -> &P {
        &self.lhs
    }

    pub fn holds(&self) -> bool {
        self.holds
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusCharacteristicEquationExhaustiveReport<P> {
    frobenius_trace: FrobeniusTrace,
    checked_points: usize,
    failed_checks: Vec<FrobeniusCharacteristicEquationCheck<P>>,
}

impl<P> FrobeniusCharacteristicEquationExhaustiveReport<P> {
    pub(crate) fn new(
        frobenius_trace: FrobeniusTrace,
        checked_points: usize,
        failed_checks: Vec<FrobeniusCharacteristicEquationCheck<P>>,
    ) -> Self {
        Self {
            frobenius_trace,
            checked_points,
            failed_checks,
        }
    }

    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    pub fn checked_points(&self) -> usize {
        self.checked_points
    }

    pub fn failed_checks(&self) -> &[FrobeniusCharacteristicEquationCheck<P>] {
        &self.failed_checks
    }

    pub fn all_hold(&self) -> bool {
        self.failed_checks.is_empty()
    }
}
