use num_bigint::BigUint;

use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;
use crate::isogenies::graphs::endomorphisms::{HorizontalEdgeReport, HorizontalEdgeStatus};

/// Compatibility status between one horizontal-edge report and one local ideal.
///
/// This status certifies only that a graph edge has already been certified
/// horizontal in a `ℓ`-volcano and that this volcano prime matches the norm of
/// a supplied prime-norm ideal. It does *not* identify a class-group action,
/// compute `E[𝔞]`, or infer which ideal should label an edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HorizontalIdealStatus {
    /// The edge is certified horizontal and the volcano prime equals the ideal norm.
    CertifiedCompatible,
    /// The edge report is not certified horizontal by crater altitude.
    EdgeNotCertifiedHorizontal,
    /// The supplied volcano prime does not equal the ideal norm.
    DegreeMismatch,
}

/// A certified compatibility witness between a horizontal edge and a prime ideal.
///
/// The witness records that the supplied edge has certified horizontal volcano evidence
/// and the supplied ideal has the same prime norm `ℓ`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HorizontalIdealWitness {
    edge: HorizontalEdgeReport,
    prime: BigUint,
    ideal: PrimeNormIdeal,
}

impl HorizontalIdealWitness {
    fn new(edge: HorizontalEdgeReport, prime: BigUint, ideal: PrimeNormIdeal) -> Self {
        Self { edge, prime, ideal }
    }

    /// Returns the certified horizontal edge report.
    pub fn edge(&self) -> &HorizontalEdgeReport {
        &self.edge
    }

    /// Returns the volcano prime `ℓ` used to certify compatibility.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the compatible prime-norm ideal.
    pub fn ideal(&self) -> &PrimeNormIdeal {
        &self.ideal
    }
}

/// Report for annotating one crater-horizontal edge with one prime-norm ideal.
///
/// The report consumes an already-computed [`HorizontalEdgeReport`], the
/// ambient volcano prime `ℓ`, and a caller-supplied [`PrimeNormIdeal`]. It
/// only checks compatibility of existing evidence:
///
/// - the edge must have status [`HorizontalEdgeStatus::CertifiedByAltitude`];
/// - `prime` must equal `ideal.norm()`.
///
/// No ideal is inferred from the edge, no orientation is exposed, and no
/// class-group action is claimed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HorizontalIdealReport {
    status: HorizontalIdealStatus,
    witness: Option<HorizontalIdealWitness>,
}

impl HorizontalIdealReport {
    /// Annotates one horizontal-edge report with one supplied prime-norm ideal.
    ///
    /// Complexity: `Θ(1)` big-integer comparison.
    pub fn from_certified_edge_and_ideal(
        edge: HorizontalEdgeReport,
        prime: BigUint,
        ideal: PrimeNormIdeal,
    ) -> Self {
        if edge.status() != HorizontalEdgeStatus::CertifiedByAltitude {
            return Self {
                status: HorizontalIdealStatus::EdgeNotCertifiedHorizontal,
                witness: None,
            };
        }

        if &prime != ideal.norm() {
            return Self {
                status: HorizontalIdealStatus::DegreeMismatch,
                witness: None,
            };
        }

        Self {
            status: HorizontalIdealStatus::CertifiedCompatible,
            witness: Some(HorizontalIdealWitness::new(edge, prime, ideal)),
        }
    }

    /// Returns the compatibility status.
    pub fn status(&self) -> HorizontalIdealStatus {
        self.status
    }

    /// Returns the witness when compatibility was certified.
    pub fn witness(&self) -> Option<&HorizontalIdealWitness> {
        self.witness.as_ref()
    }
}
