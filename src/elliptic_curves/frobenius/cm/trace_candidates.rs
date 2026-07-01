use core::fmt;

use num_bigint::BigUint;
use num_traits::{Signed, Zero};

use crate::elliptic_curves::endomorphisms::quadratic_orders::QuadraticDiscriminant;
use crate::numerics::cornacchia::{CornacchiaError, cornacchia_candidate_solutions};
use crate::numerics::quadratic_forms::{DiagonalBinaryQuadraticForm, QuadraticFormError};

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

/// Failure modes for CM trace-candidate generation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CmTraceCandidateError {
    /// The CM discriminant must be negative.
    NonNegativeDiscriminant,
    /// The prime-like input `p` must be positive.
    ZeroPrime,
    /// The underlying Cornacchia candidate route failed.
    Cornacchia(CornacchiaError),
    /// The auxiliary primitive representation route for `p` failed.
    QuadraticForm(QuadraticFormError),
}

impl From<CornacchiaError> for CmTraceCandidateError {
    fn from(error: CornacchiaError) -> Self {
        Self::Cornacchia(error)
    }
}

impl From<QuadraticFormError> for CmTraceCandidateError {
    fn from(error: QuadraticFormError) -> Self {
        Self::QuadraticForm(error)
    }
}

impl fmt::Display for CmTraceCandidateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonNegativeDiscriminant => write!(
                formatter,
                "CM trace candidates require a negative discriminant D"
            ),
            Self::ZeroPrime => write!(
                formatter,
                "CM trace candidates require a positive prime-like input p"
            ),
            Self::Cornacchia(error) => {
                write!(formatter, "Cornacchia candidate generation failed: {error}")
            }
            Self::QuadraticForm(error) => {
                write!(formatter, "quadratic-form representation failed: {error}")
            }
        }
    }
}

impl std::error::Error for CmTraceCandidateError {}

/// Returns candidate absolute traces `|t|` from `4p = t² + |D|v²`.
///
/// The input `discriminant` is interpreted as the CM discriminant `D < 0`.
/// This helper forms the diagonal equation `x² + |D|y² = 4p` and returns the
/// represented `x` values as candidates for `|t|`.
///
/// The implementation deliberately combines two sources:
///
/// - direct Cornacchia candidates for `4p = x² + |D|y²`, and
/// - primitive representations `p = u² + |D|w²`, lifted to
///   `4p = (2u)² + |D|(2w)²`.
///
/// This is an arithmetic Frobenius-side helper, not a curve-CM certificate:
/// it does not prove that any concrete curve has CM by `D`, and it intentionally
/// does not determine the sign of the Frobenius trace.
///
/// Complexity: dominated by candidate generation for `4p = x² + |D|y²` and
/// primitive representation of `p` by `x² + |D|y²`.
pub fn cm_absolute_trace_candidates(
    discriminant: &QuadraticDiscriminant,
    p: &BigUint,
) -> Result<Vec<BigUint>, CmTraceCandidateError> {
    Ok(cm_trace_candidates(discriminant, p)?
        .into_iter()
        .map(|candidate| candidate.absolute_trace)
        .collect())
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

    let d = absolute_discriminant(discriminant);
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

fn absolute_discriminant(discriminant: &QuadraticDiscriminant) -> BigUint {
    (-discriminant.value())
        .to_biguint()
        .expect("negating a negative discriminant should produce a positive integer")
}
