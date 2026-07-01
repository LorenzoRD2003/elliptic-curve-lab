use num_bigint::BigUint;
use num_traits::{Signed, Zero};

use crate::elliptic_curves::{
    endomorphisms::quadratic_orders::QuadraticDiscriminant, frobenius::cm::CmTraceCandidateError,
};
use crate::numerics::{
    cornacchia::cornacchia_candidate_solutions, quadratic_forms::DiagonalBinaryQuadraticForm,
};

/// One arithmetic CM candidate satisfying `4p = |t|² + |D|v²`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct CmTraceCandidate {
    absolute_trace: BigUint,
    cm_multiplier: BigUint,
}

impl CmTraceCandidate {
    fn new(absolute_trace: BigUint, cm_multiplier: BigUint) -> Self {
        Self {
            absolute_trace,
            cm_multiplier,
        }
    }

    /// Returns the candidate absolute Frobenius trace `|t|`.
    pub fn absolute_trace(&self) -> &BigUint {
        &self.absolute_trace
    }

    /// Returns the auxiliary multiplier `v` in `4p = |t|² + |D|v²`.
    pub fn cm_multiplier(&self) -> &BigUint {
        &self.cm_multiplier
    }
}

/// Returns candidate absolute traces with witnesses from `4p = t² + |D|v²`.
///
/// The input `discriminant` is interpreted as the CM discriminant `D < 0`.
/// This helper forms the diagonal equation `x² + |D|y² = 4p` and returns
/// `(x, y)` as [`CmTraceCandidate`] values, where `x = |t|` and `y = v`.
///
/// The implementation combines two sources:
///
/// - direct Cornacchia candidates for `4p = x² + |D|y²`, and
/// - primitive representations `p = u² + |D|w²`, lifted to
///   `4p = (2u)² + |D|(2w)²`.
///
/// This does *not* prove that any concrete curve has CM by `D`, and it
/// does *not* determine the sign of the Frobenius trace.
///
/// Complexity: `Θ(C_4p + C_p + k log k)`, where `C_4p` is the cost of
/// Cornacchia candidate generation for `4p = x² + |D|y²`, `C_p` is the cost of
/// primitive representation of `p` by `x² + |D|y²`, and `k` is the number of
/// candidates before deduplication.
pub fn cm_absolute_trace_candidates(
    discriminant: &QuadraticDiscriminant,
    p: &BigUint,
) -> Result<Vec<CmTraceCandidate>, CmTraceCandidateError> {
    cm_trace_candidates(discriminant, p)
}

fn cm_trace_candidates(
    discriminant: &QuadraticDiscriminant,
    p: &BigUint,
) -> Result<Vec<CmTraceCandidate>, CmTraceCandidateError> {
    if !discriminant.value().is_negative() {
        return Err(CmTraceCandidateError::NonNegativeDiscriminant);
    }
    if p.is_zero() {
        return Err(CmTraceCandidateError::ZeroPrime);
    }

    let d = discriminant.value().magnitude().clone();
    let m = BigUint::from(4u8) * p;
    let mut candidates = cornacchia_candidate_solutions(&d, &m)?
        .into_iter()
        .map(|solution| CmTraceCandidate::new(solution.x().clone(), solution.y().clone()))
        .collect::<Vec<_>>();

    let form = DiagonalBinaryQuadraticForm::new(d)?;
    candidates.extend(
        form.primitive_representations(p)?
            .into_iter()
            .map(|representation| {
                CmTraceCandidate::new(
                    BigUint::from(2u8) * representation.x(),
                    BigUint::from(2u8) * representation.y(),
                )
            }),
    );

    candidates.sort();
    candidates.dedup();
    Ok(candidates)
}
