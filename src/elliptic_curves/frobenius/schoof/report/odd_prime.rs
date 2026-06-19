use num_bigint::BigUint;

use crate::elliptic_curves::short_weierstrass::{
    ReducedEndomorphism, ReducedEndomorphismAdditiveResult,
};
use crate::fields::traits::FiniteField;
use crate::numerics::chinese_remainder::Congruence;
use crate::polynomials::DensePolynomial;

/// One tested candidate `c mod ℓ` in the odd-prime Schoof trace loop.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofTraceModOddPrimeCandidateReport<F: FiniteField> {
    candidate_trace_mod_ell: usize,
    result: ReducedEndomorphismAdditiveResult<F>,
}

impl<F: FiniteField> SchoofTraceModOddPrimeCandidateReport<F> {
    pub(crate) fn new(
        candidate_trace_mod_ell: usize,
        result: ReducedEndomorphismAdditiveResult<F>,
    ) -> Self {
        Self {
            candidate_trace_mod_ell,
            result,
        }
    }

    /// Returns the tested candidate residue `c mod ℓ`.
    pub fn candidate_trace_mod_ell(&self) -> usize {
        self.candidate_trace_mod_ell
    }

    /// Returns the reduced characteristic-equation result for this candidate.
    pub fn result(&self) -> &ReducedEndomorphismAdditiveResult<F> {
        &self.result
    }
}

/// Final outcome of the current odd-prime Schoof trace driver.
#[derive(Clone, Debug, PartialEq)]
pub enum SchoofTraceModOddPrimeOutcome<F: FiniteField> {
    /// The driver found the unique candidate `c mod ℓ` for which
    /// `π^2 - [c]π + [q] id = 0` in the current reduced arithmetic.
    TraceFound { trace_mod_ell: usize },
    /// The driver stopped early because one denominator was not invertible
    /// modulo `ψ_ℓ(x)`.
    NonUnitDenominator {
        candidate_trace_mod_ell: usize,
        witness_gcd: DensePolynomial<F>,
    },
    /// The driver tested every residue class modulo `ℓ` and did not find a
    /// zero candidate before refinement-by-factor was attempted.
    ExhaustedCandidates,
}

/// Report for the odd-prime `ℓ` step of Schoof's trace computation.
#[derive(Clone, Debug, PartialEq)]
pub struct SchoofTraceModOddPrimeReport<F: FiniteField> {
    field_order: u128,
    odd_prime: usize,
    division_polynomial: DensePolynomial<F>,
    frobenius: ReducedEndomorphism<F>,
    frobenius_squared: ReducedEndomorphism<F>,
    candidate_reports: Vec<SchoofTraceModOddPrimeCandidateReport<F>>,
    outcome: SchoofTraceModOddPrimeOutcome<F>,
}

impl<F: FiniteField> SchoofTraceModOddPrimeReport<F> {
    pub(crate) fn new(
        field_order: u128,
        odd_prime: usize,
        division_polynomial: DensePolynomial<F>,
        frobenius: ReducedEndomorphism<F>,
        frobenius_squared: ReducedEndomorphism<F>,
        candidate_reports: Vec<SchoofTraceModOddPrimeCandidateReport<F>>,
        outcome: SchoofTraceModOddPrimeOutcome<F>,
    ) -> Self {
        Self {
            field_order,
            odd_prime,
            division_polynomial,
            frobenius,
            frobenius_squared,
            candidate_reports,
            outcome,
        }
    }

    /// Returns the finite field order `q`.
    pub fn field_order(&self) -> u128 {
        self.field_order
    }

    /// Returns the odd prime `ℓ`.
    pub fn odd_prime(&self) -> usize {
        self.odd_prime
    }

    /// Returns the odd division polynomial `ψ_ℓ(x)` used as quotient modulus.
    pub fn division_polynomial(&self) -> &DensePolynomial<F> {
        &self.division_polynomial
    }

    /// Returns the reduced Frobenius endomorphism `π`.
    pub fn frobenius(&self) -> &ReducedEndomorphism<F> {
        &self.frobenius
    }

    /// Returns the reduced composed endomorphism `π²`.
    pub fn frobenius_squared(&self) -> &ReducedEndomorphism<F> {
        &self.frobenius_squared
    }

    /// Returns the tested candidate reports, in loop order.
    pub fn candidate_reports(&self) -> &[SchoofTraceModOddPrimeCandidateReport<F>] {
        &self.candidate_reports
    }

    /// Returns the final driver outcome.
    pub fn outcome(&self) -> &SchoofTraceModOddPrimeOutcome<F> {
        &self.outcome
    }

    /// Returns `tr(π_q) mod ℓ` when the current driver already found it
    /// without factor refinement.
    pub fn trace_mod_odd_prime(&self) -> Option<usize> {
        match self.outcome {
            SchoofTraceModOddPrimeOutcome::TraceFound { trace_mod_ell } => Some(trace_mod_ell),
            _ => None,
        }
    }

    /// Returns the congruence `t ≡ c (mod ℓ)` when this report already found
    /// one trace residue `c mod ℓ` without factor refinement.
    pub fn trace_congruence(&self) -> Option<Congruence> {
        self.trace_mod_odd_prime().map(|trace_mod_ell| {
            Congruence::new(
                BigUint::from(trace_mod_ell),
                BigUint::from(self.odd_prime as u64),
            )
            .expect("an odd prime modulus should define a valid congruence")
        })
    }
}
