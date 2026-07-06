use std::collections::BTreeSet;

use num_bigint::BigUint;
use num_traits::One;

use crate::elliptic_curves::endomorphisms::{
    candidate_sets::EndomorphismRingCandidateSet, quadratic_orders::ImaginaryQuadraticOrder,
};
use crate::isogenies::graphs::{
    IsogenyGraphNodeId,
    endomorphisms::{EndomorphismRingLevelRecoveryError, LocalEndomorphismRingLevelReport},
};
use crate::numerics::distinct_prime_factors;

/// Global endomorphism-ring recovery assembled from local volcano reports.
///
/// If `Δ_π = v²D_K`, then every candidate endomorphism ring has the form
/// `O_u` with `u | v`. For each prime factor `ℓᵉ || v`, a local report recovers
/// `d = v_ℓ(u)` from the volcano-floor identity `d = e - δ`.
///
/// This report combines those local exponents. It constructs the recovered
/// order only when the supplied reports cover every prime divisor of `v`;
/// otherwise it remains a partial report and exposes the missing primes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EndomorphismRingLevelRecoveryReport {
    candidate_set: EndomorphismRingCandidateSet,
    local_reports: Vec<LocalEndomorphismRingLevelReport>,
    missing_primes: Vec<BigUint>,
    recovered_conductor: Option<BigUint>,
    recovered_order: Option<ImaginaryQuadraticOrder>,
    node_id: Option<IsogenyGraphNodeId>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ValidatedLocalReports {
    node_id: Option<IsogenyGraphNodeId>,
    seen_primes: BTreeSet<BigUint>,
}

impl EndomorphismRingLevelRecoveryReport {
    /// Builds a global recovery report from local `ℓ`-level reports.
    ///
    /// The method validates that the local reports:
    ///
    /// 1. all refer to the same graph node,
    /// 2. use distinct local primes,
    /// 3. only use primes dividing the Frobenius conductor `v`,
    /// 4. agree with the candidate set on `e = v_ℓ(v)`.
    ///
    /// If every prime divisor of `v` is present, the report reconstructs
    /// `u = ∏ℓ^{d_ℓ}` and the order `O_u`. If at least one prime is missing,
    /// [`Self::recovered_conductor`] and [`Self::recovered_order`] return
    /// `None`.
    ///
    /// Complexity: dominated by factoring the Frobenius conductor `v` with
    /// `num-prime`, plus linear validation in the number of local reports.
    pub fn from_local_reports(
        candidate_set: EndomorphismRingCandidateSet,
        mut local_reports: Vec<LocalEndomorphismRingLevelReport>,
    ) -> Result<Self, EndomorphismRingLevelRecoveryError> {
        local_reports.sort_by(|left, right| left.prime().cmp(right.prime()));

        let required_primes = distinct_prime_factors(candidate_set.frobenius_conductor());
        let validated =
            Self::validate_local_reports(&candidate_set, &required_primes, &local_reports)?;
        let missing_primes = required_primes
            .iter()
            .filter(|prime| !validated.seen_primes.contains(*prime))
            .cloned()
            .collect::<Vec<_>>();

        let (recovered_conductor, recovered_order) = if missing_primes.is_empty() {
            let conductor = Self::recovered_conductor_from_local_reports(&local_reports);
            let order = candidate_set
                .candidate_orders()
                .iter()
                .find(|order| order.conductor() == &conductor)
                .expect("recovered conductor is assembled from divisors of v")
                .clone();
            (Some(conductor), Some(order))
        } else {
            (None, None)
        };

        Ok(Self {
            candidate_set,
            local_reports,
            missing_primes,
            recovered_conductor,
            recovered_order,
            node_id: validated.node_id,
        })
    }

    /// Returns the node shared by the local reports, when at least one local
    /// report was supplied.
    pub fn node_id(&self) -> Option<IsogenyGraphNodeId> {
        self.node_id
    }

    /// Returns the Frobenius-compatible candidate set `C₀`.
    pub fn candidate_set(&self) -> &EndomorphismRingCandidateSet {
        &self.candidate_set
    }

    /// Returns the local reports, sorted increasingly by prime `ℓ`.
    pub fn local_reports(&self) -> &[LocalEndomorphismRingLevelReport] {
        &self.local_reports
    }

    /// Returns the recovered conductor `u` when every prime `ℓ | v` is covered.
    pub fn recovered_conductor(&self) -> Option<&BigUint> {
        self.recovered_conductor.as_ref()
    }

    /// Returns the recovered order `O_u` when every prime `ℓ | v` is covered.
    pub fn recovered_order(&self) -> Option<&ImaginaryQuadraticOrder> {
        self.recovered_order.as_ref()
    }

    /// Returns the primes `ℓ | v` still missing local recovery evidence.
    pub fn missing_primes(&self) -> &[BigUint] {
        &self.missing_primes
    }

    /// Returns whether the local reports cover every prime divisor of `v`.
    pub fn is_complete(&self) -> bool {
        self.missing_primes.is_empty()
    }

    /// Validates that the supplied local reports can be assembled against this
    /// candidate set.
    ///
    /// This helper asserts the following invariants: each local report must describe
    /// the same node, a distinct prime `ℓ` dividing the Frobenius conductor `v`,
    /// and the exponent `e_ℓ = v_ℓ(v)` that the candidate set derives from `Δ_π = v²D_K`.
    ///
    /// The returned [`ValidatedLocalReports`] records the shared node id and
    /// the set of observed primes, which the caller then compares with the full
    /// set of prime divisors of `v` to decide if the global recovery is complete.
    fn validate_local_reports(
        candidate_set: &EndomorphismRingCandidateSet,
        required_primes: &[BigUint],
        local_reports: &[LocalEndomorphismRingLevelReport],
    ) -> Result<ValidatedLocalReports, EndomorphismRingLevelRecoveryError> {
        let mut node_id = None;
        let mut seen_primes = BTreeSet::new();

        for report in local_reports {
            match node_id {
                Some(expected_node_id) if expected_node_id != report.node_id() => {
                    return Err(EndomorphismRingLevelRecoveryError::MixedNodeReports {
                        expected_node_id,
                        found_node_id: report.node_id(),
                    });
                }
                None => node_id = Some(report.node_id()),
                Some(_) => {}
            }

            let prime = report.prime();
            if !seen_primes.insert(prime.clone()) {
                return Err(EndomorphismRingLevelRecoveryError::DuplicateLocalPrime {
                    prime: prime.clone(),
                });
            }

            if !required_primes.contains(prime) {
                return Err(
                    EndomorphismRingLevelRecoveryError::LocalPrimeNotInFrobeniusConductor {
                        prime: prime.clone(),
                    },
                );
            }

            let expected = candidate_set
                .local_view_at(prime)?
                .frobenius_conductor_valuation();
            if report.frobenius_conductor_valuation() != expected {
                return Err(
                    EndomorphismRingLevelRecoveryError::InconsistentLocalConductorValuation {
                        prime: prime.clone(),
                        report_frobenius_conductor_valuation: report
                            .frobenius_conductor_valuation(),
                        expected_frobenius_conductor_valuation: expected,
                    },
                );
            }
        }

        Ok(ValidatedLocalReports {
            node_id,
            seen_primes,
        })
    }

    fn recovered_conductor_from_local_reports(
        local_reports: &[LocalEndomorphismRingLevelReport],
    ) -> BigUint {
        local_reports
            .iter()
            .fold(BigUint::one(), |conductor, report| {
                conductor * report.prime().pow(report.recovered_conductor_valuation())
            })
    }
}
