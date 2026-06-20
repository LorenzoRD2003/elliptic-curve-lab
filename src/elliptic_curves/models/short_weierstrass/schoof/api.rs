use num_bigint::BigUint;
use num_prime::nt_funcs::is_prime;

use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::{
        HasseInterval,
        schoof::{
            SchoofGroupOrderReport, SchoofTraceCrtOutcome, SchoofTraceCrtReport,
            SchoofTraceMod2Report, SchoofTraceModOddPrimeCandidateReport,
            SchoofTraceModOddPrimeOutcome, SchoofTraceModOddPrimeReport,
            finalize_schoof_group_order_report,
        },
    },
    short_weierstrass::{
        division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
        schoof::{ReducedCurveQuotient, ReducedEndomorphismAdditiveResult},
    },
};
use crate::fields::{finite_field_descriptor::FiniteFieldDescriptor, traits::FiniteField};
use crate::numerics::chinese_remainder::{ChineseRemainderSolution, combine_coprime_congruences};
use crate::polynomials::DensePolynomial;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Computes `tr(π_q) mod 2` through the first step of Schoof's algorithm.
    ///
    /// For odd `q`, this determines whether `E(F_q)` has rational `2`-torsion
    /// by checking whether the short-Weierstrass cubic has a root in `F_q`.
    ///
    /// The implementation computes `x^q mod f(x)` first and then forms
    /// `gcd(f(x), (x^q mod f(x)) - x)`, matching the quotient-ring optimization
    /// already used elsewhere in the crate.
    pub fn schoof_trace_mod_2(&self) -> SchoofTraceMod2Report<F> {
        let cubic = self.to_cubic();
        let x = DensePolynomial::new(vec![F::zero(), F::one()]);
        let field_order = F::order().expect("finite field order should fit in u128");
        let x_q_mod_cubic = DensePolynomial::pow_mod(&x, field_order, &cubic)
            .expect("short-Weierstrass cubic is a non-zero modulus");
        let gcd = cubic.gcd(&x_q_mod_cubic.sub(&x));

        SchoofTraceMod2Report::new(field_order, cubic, x_q_mod_cubic, gcd)
    }

    /// Runs the current Schoof prime steps and combines every resolved trace
    /// congruence by the Chinese remainder theorem.
    ///
    /// The current implementation always starts with `ℓ = 2`, then processes
    /// the requested odd primes in order. Each odd-prime step is the current
    /// non-refining driver: if one such step stops at a non-unit denominator,
    /// this CRT stage also stops and returns the partial CRT solution
    /// accumulated strictly before that prime.
    ///
    /// Complexity: the sum of the current `ℓ = 2` and odd-prime Schoof step
    /// costs, plus `Θ(r)` CRT combinations and `Θ(r^2)` pairwise-coprimality
    /// checks for `r` resolved congruences.
    pub fn schoof_trace_crt(
        &self,
        odd_primes: &[usize],
    ) -> Result<SchoofTraceCrtReport<F>, DivisionPolynomialError> {
        let mut state = SchoofCrtState::new(self.schoof_trace_mod_2());
        for &odd_prime in odd_primes {
            state = match self.extend_schoof_crt_state(state, odd_prime)? {
                SchoofCrtExtension::Continue(state) => state,
                SchoofCrtExtension::Skipped {
                    next_state,
                    skipped_prime,
                } => return Ok(next_state.blocked_on_odd_prime(skipped_prime)),
            };
        }
        Ok(state.finish())
    }

    /// Runs Schoof's CRT accumulation with the natural stopping condition:
    /// keep adding odd primes `ℓ = 3, 5, 7, ...` until the combined CRT
    /// modulus exceeds the Hasse trace diameter bound `2⌊2√q⌋`.
    ///
    /// At that point, any residue class modulo the accumulated modulus can
    /// contain at most one trace inside `[-⌊2√q⌋, ⌊2√q⌋]`, so the final
    /// Hasse-resolution stage can no longer be ambiguous.
    ///
    /// If one odd-prime step blocks on a non-unit denominator before that
    /// threshold is reached, the current automatic route records that prime's
    /// report but skips its congruence contribution and continues with the next
    /// odd prime instead of aborting.
    ///
    /// Complexity: the sum of the current Schoof prime-step costs for the
    /// attempted primes, plus `Θ(r)` CRT combinations for `r` resolved
    /// congruences before the route either blocks or crosses the Hasse
    /// uniqueness threshold.
    pub fn schoof_trace_crt_until_hasse_uniqueness(
        &self,
    ) -> Result<SchoofTraceCrtReport<F>, DivisionPolynomialError> {
        let mut state = SchoofCrtState::new(self.schoof_trace_mod_2());
        let uniqueness_threshold = schoof_trace_uniqueness_threshold(state.field_order)?;
        if state.partial_solution.modulus() > &uniqueness_threshold {
            return Ok(state.finish());
        }

        let mut odd_prime = 3usize;
        loop {
            state = match self.extend_schoof_crt_state(state, odd_prime)? {
                SchoofCrtExtension::Continue(next_state) => next_state,
                SchoofCrtExtension::Skipped {
                    next_state,
                    skipped_prime: _,
                } => next_state,
            };
            if state.partial_solution.modulus() > &uniqueness_threshold {
                return Ok(state.finish());
            }
            odd_prime = next_schoof_odd_prime_after(odd_prime, F::characteristic());
        }
    }

    /// Runs the natural end-to-end Schoof route that keeps adding odd primes
    /// until the accumulated CRT modulus is large enough for Hasse's bound to
    /// force one unique Frobenius trace.
    ///
    /// This is the intended public surface for the general finite-field
    /// group-order route whose `Auto` policy is Schoof. The manual odd-prime list remains
    /// available only on the trace/CRT side as an educational inspection
    /// surface.
    ///
    /// Complexity: the cost of
    /// [`Self::schoof_trace_crt_until_hasse_uniqueness`] plus `Θ(1)` exact
    /// integer arithmetic for the final Hasse-resolution step.
    pub fn schoof_group_order(&self) -> Result<SchoofGroupOrderReport<F>, DivisionPolynomialError> {
        let crt_report = self.schoof_trace_crt_until_hasse_uniqueness()?;
        let base_field = FiniteFieldDescriptor::new(F::characteristic(), F::extension_degree())
            .map_err(|_| CurveError::InvalidFrobeniusBaseField {
                characteristic: F::characteristic(),
                extension_degree: F::extension_degree().get(),
            })?;
        finalize_schoof_group_order_report(base_field, crt_report).map_err(Into::into)
    }

    /// Runs the current odd-prime `ℓ` step of Schoof's algorithm *without*
    /// refinement-by-factor fallback.
    ///
    /// - validates that `ℓ` is an odd prime different from `char(F_q)`
    /// - builds the odd division polynomial `ψ_ℓ(x)`
    /// - works in the reduced quotient `F[x, y] / (y^2 - f(x), ψ_ℓ(x))`
    /// - constructs `π`, `π²`, and tests `π² - [c]π + [q] id` for `c = 0, ..., ℓ - 1`
    ///
    /// If one candidate yields the additive zero endomorphism, this returns
    /// that trace residue modulo `ℓ`. If a denominator becomes noninvertible
    /// modulo `ψ_ℓ(x)`, the current implementation stops and records the gcd
    /// witness rather than refining by factors.
    ///
    /// Complexity: `Θ(ℓ (log ℓ + log q + m) m^2)` field operations, where
    /// `m = deg ψ_ℓ`, dominated by testing up to `ℓ` characteristic-equation
    /// candidates in the reduced additive arithmetic.
    pub fn schoof_trace_mod_odd_prime(
        &self,
        odd_prime: usize,
    ) -> Result<SchoofTraceModOddPrimeReport<F>, DivisionPolynomialError> {
        self.ensure_valid_schoof_odd_prime(odd_prime)?;

        let division_polynomial = match self.division_polynomial(odd_prime)? {
            DivisionPolynomialForm::InX(polynomial) => polynomial,
            DivisionPolynomialForm::YTimes(_) => {
                return Err(CurveError::InvalidSchoofOddPrime {
                    odd_prime,
                    characteristic: F::characteristic(),
                }
                .into());
            }
        };

        let field_order = F::order().expect("finite field order should fit in u128");
        let quotient = ReducedCurveQuotient::new(self.clone(), division_polynomial.clone())?;
        let frobenius = self.reduced_frobenius_endomorphism(&quotient);
        let frobenius_squared = frobenius.compose(&quotient, &frobenius);
        let q_term = self.scalar_multiple_of_reduced_identity_endomorphism_on_odd_torsion(
            &quotient, odd_prime, field_order,
        );

        let mut candidate_reports = Vec::with_capacity(odd_prime);

        for candidate_trace_mod_ell in 0..odd_prime {
            let result = self.reduced_characteristic_equation_candidate(
                &quotient,
                odd_prime,
                &frobenius,
                &frobenius_squared,
                &q_term,
                candidate_trace_mod_ell as u128,
            );
            candidate_reports.push(SchoofTraceModOddPrimeCandidateReport::new(
                candidate_trace_mod_ell,
                result.clone(),
            ));

            match result {
                ReducedEndomorphismAdditiveResult::Zero => {
                    return Ok(SchoofTraceModOddPrimeReport::new(
                        field_order,
                        odd_prime,
                        division_polynomial,
                        frobenius.clone(),
                        frobenius_squared.clone(),
                        candidate_reports,
                        SchoofTraceModOddPrimeOutcome::TraceFound {
                            trace_mod_ell: candidate_trace_mod_ell,
                        },
                    ));
                }
                ReducedEndomorphismAdditiveResult::NonUnitDenominator { witness_gcd } => {
                    // TODO: refine the active modulus `ψ_ℓ(x)` by this nontrivial gcd
                    // witness and continue the odd-prime Schoof step recursively on the
                    // resulting factors, instead of stopping at the first non-unit branch.
                    return Ok(SchoofTraceModOddPrimeReport::new(
                        field_order,
                        odd_prime,
                        division_polynomial,
                        frobenius.clone(),
                        frobenius_squared.clone(),
                        candidate_reports,
                        SchoofTraceModOddPrimeOutcome::NonUnitDenominator {
                            candidate_trace_mod_ell,
                            witness_gcd,
                        },
                    ));
                }
                ReducedEndomorphismAdditiveResult::Value(_) => {}
            }
        }

        Ok(SchoofTraceModOddPrimeReport::new(
            field_order,
            odd_prime,
            division_polynomial,
            frobenius,
            frobenius_squared,
            candidate_reports,
            SchoofTraceModOddPrimeOutcome::ExhaustedCandidates,
        ))
    }

    fn ensure_valid_schoof_odd_prime(
        &self,
        odd_prime: usize,
    ) -> Result<(), DivisionPolynomialError> {
        if odd_prime == 2
            || !is_prime(&(odd_prime as u64), None).probably()
            || odd_prime as u64 == F::characteristic()
        {
            Err(CurveError::InvalidSchoofOddPrime {
                odd_prime,
                characteristic: F::characteristic(),
            }
            .into())
        } else {
            Ok(())
        }
    }

    fn extend_schoof_crt_state(
        &self,
        mut state: SchoofCrtState<F>,
        odd_prime: usize,
    ) -> Result<SchoofCrtExtension<F>, DivisionPolynomialError> {
        let report = self.schoof_trace_mod_odd_prime(odd_prime)?;
        let trace_congruence = report.trace_congruence();
        state.odd_prime_reports.push(report);

        let Some(congruence) = trace_congruence else {
            return Ok(SchoofCrtExtension::Skipped {
                next_state: state,
                skipped_prime: odd_prime,
            });
        };

        state.partial_solution = combine_coprime_congruences(&state.partial_solution, &congruence)
            .expect("distinct Schoof primes should stay coprime in CRT recombination");
        state.resolved_congruences.push(congruence);
        Ok(SchoofCrtExtension::Continue(state))
    }
}

struct SchoofCrtState<F: FiniteField> {
    field_order: u128,
    mod_2_report: SchoofTraceMod2Report<F>,
    odd_prime_reports: Vec<SchoofTraceModOddPrimeReport<F>>,
    resolved_congruences: Vec<crate::numerics::chinese_remainder::Congruence>,
    partial_solution: ChineseRemainderSolution,
}

impl<F: FiniteField> SchoofCrtState<F> {
    fn new(mod_2_report: SchoofTraceMod2Report<F>) -> Self {
        let field_order = mod_2_report.field_order();
        let mod_2_congruence = mod_2_report.trace_congruence();
        let partial_solution = ChineseRemainderSolution::new(
            mod_2_congruence.residue().clone(),
            mod_2_congruence.modulus().clone(),
        );

        Self {
            field_order,
            mod_2_report,
            odd_prime_reports: Vec::new(),
            resolved_congruences: vec![mod_2_congruence],
            partial_solution,
        }
    }

    fn blocked_on_odd_prime(self, blocked_prime: usize) -> SchoofTraceCrtReport<F> {
        SchoofTraceCrtReport::new(
            self.field_order,
            self.mod_2_report,
            self.odd_prime_reports,
            self.resolved_congruences,
            SchoofTraceCrtOutcome::BlockedOnOddPrime {
                blocked_prime,
                partial_solution: self.partial_solution,
            },
        )
    }

    fn finish(self) -> SchoofTraceCrtReport<F> {
        SchoofTraceCrtReport::new(
            self.field_order,
            self.mod_2_report,
            self.odd_prime_reports,
            self.resolved_congruences,
            SchoofTraceCrtOutcome::Combined {
                solution: self.partial_solution,
            },
        )
    }
}

enum SchoofCrtExtension<F: FiniteField> {
    Continue(SchoofCrtState<F>),
    Skipped {
        next_state: SchoofCrtState<F>,
        skipped_prime: usize,
    },
}

fn schoof_trace_uniqueness_threshold(field_order: u128) -> Result<BigUint, CurveError> {
    let interval = HasseInterval::for_q(field_order)?;
    Ok(BigUint::from(interval.trace_bound()) * BigUint::from(2u8))
}

fn next_schoof_odd_prime_after(previous_odd_prime: usize, characteristic: u64) -> usize {
    let mut candidate = previous_odd_prime
        .checked_add(2)
        .expect("educational Schoof prime search should not exhaust usize");
    loop {
        if is_prime(&(candidate as u64), None).probably() && candidate as u64 != characteristic {
            return candidate;
        }
        candidate = candidate
            .checked_add(2)
            .expect("educational Schoof prime search should not exhaust usize");
    }
}
