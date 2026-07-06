use std::fmt;

use num_bigint::BigUint;

use crate::isogenies::graphs::{
    IsogenyGraphError, IsogenyGraphNodeId, endomorphisms::VolcanoSearchError,
};
use crate::numerics::PositivePrimeError;

/// Errors produced while recovering one local endomorphism-ring level.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EndomorphismRingLevelRecoveryError {
    /// The chosen local parameter is not a positive prime `ℓ`.
    InvalidLocalPrime(PositivePrimeError),
    /// Frobenius-compatible candidate-order derivation failed for the node.
    CandidateDerivation(IsogenyGraphError),
    /// Floor-distance certification failed for the node in the `ℓ`-graph.
    FloorSearch(VolcanoSearchError),
    /// The certified floor distance is larger than `v_ℓ(v)`.
    ///
    /// For an ordinary `ℓ`-volcano compatible with `Δ_π = v²D_K`, the formula
    /// `d = e - δ` requires `δ <= e`, where `e = v_ℓ(v)`.
    DistanceExceedsFrobeniusConductorValuation {
        /// The node being recovered.
        node_id: IsogenyGraphNodeId,
        /// The chosen local prime `ℓ`.
        prime: BigUint,
        /// The certified distance `δ` from the node to the floor.
        distance_to_floor: usize,
        /// The Frobenius conductor valuation `e = v_ℓ(v)`.
        frobenius_conductor_valuation: u32,
    },
    /// A global recovery report received two local reports for the same prime.
    DuplicateLocalPrime {
        /// The duplicated local prime `ℓ`.
        prime: BigUint,
    },
    /// A local report was supplied for a prime not dividing `v`.
    LocalPrimeNotInFrobeniusConductor {
        /// The local prime from the report.
        prime: BigUint,
    },
    /// A local report's `e = v_ℓ(v)` disagrees with the global candidate set.
    InconsistentLocalConductorValuation {
        /// The chosen local prime `ℓ`.
        prime: BigUint,
        /// The exponent recorded in the local report.
        report_frobenius_conductor_valuation: u32,
        /// The exponent derived from the global candidate set.
        expected_frobenius_conductor_valuation: u32,
    },
    /// Local reports from different nodes were combined into one global run.
    MixedNodeReports {
        /// Node id established by the first local report.
        expected_node_id: IsogenyGraphNodeId,
        /// A later local report's node id.
        found_node_id: IsogenyGraphNodeId,
    },
}

impl fmt::Display for EndomorphismRingLevelRecoveryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLocalPrime(error) => {
                write!(formatter, "invalid local recovery prime: {error}")
            }
            Self::CandidateDerivation(error) => write!(
                formatter,
                "could not derive Frobenius-compatible endomorphism-ring candidates: {error}"
            ),
            Self::FloorSearch(error) => write!(
                formatter,
                "could not certify the volcano floor distance needed for local recovery: {error}"
            ),
            Self::DistanceExceedsFrobeniusConductorValuation {
                node_id,
                prime,
                distance_to_floor,
                frobenius_conductor_valuation,
            } => write!(
                formatter,
                "node {:?} has certified ℓ-floor distance {distance_to_floor} at ℓ = {prime}, larger than v_ℓ(v) = {frobenius_conductor_valuation}",
                node_id
            ),
            Self::DuplicateLocalPrime { prime } => write!(
                formatter,
                "global endomorphism-ring recovery received more than one local report for ℓ = {prime}"
            ),
            Self::LocalPrimeNotInFrobeniusConductor { prime } => write!(
                formatter,
                "local report at ℓ = {prime} cannot contribute to global recovery because ℓ does not divide the Frobenius conductor v"
            ),
            Self::InconsistentLocalConductorValuation {
                prime,
                report_frobenius_conductor_valuation,
                expected_frobenius_conductor_valuation,
            } => write!(
                formatter,
                "local report at ℓ = {prime} records v_ℓ(v) = {report_frobenius_conductor_valuation}, but the global candidate set has v_ℓ(v) = {expected_frobenius_conductor_valuation}"
            ),
            Self::MixedNodeReports {
                expected_node_id,
                found_node_id,
            } => write!(
                formatter,
                "global endomorphism-ring recovery received reports from different nodes: expected {:?}, found {:?}",
                expected_node_id, found_node_id
            ),
        }
    }
}

impl std::error::Error for EndomorphismRingLevelRecoveryError {}

impl From<PositivePrimeError> for EndomorphismRingLevelRecoveryError {
    fn from(error: PositivePrimeError) -> Self {
        Self::InvalidLocalPrime(error)
    }
}

impl From<IsogenyGraphError> for EndomorphismRingLevelRecoveryError {
    fn from(error: IsogenyGraphError) -> Self {
        Self::CandidateDerivation(error)
    }
}

impl From<VolcanoSearchError> for EndomorphismRingLevelRecoveryError {
    fn from(error: VolcanoSearchError) -> Self {
        Self::FloorSearch(error)
    }
}
