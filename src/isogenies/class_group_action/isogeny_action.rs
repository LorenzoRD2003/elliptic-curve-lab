use core::fmt;

use num_bigint::{BigInt, BigUint};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::BinaryQuadraticForm, quadratic_ideals::PrimeNormIdeal,
};
use crate::isogenies::{
    class_group_action::{CraterDirectionCertification, OrientedCraterPowerActionError},
    graphs::IsogenyGraphNodeId,
};

/// Failure modes for executing an algebraic class-group action plan on oriented craters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ClassGroupIsogenyActionError {
    /// No oriented local crater witness matched the requested plan factor.
    MissingLocalWitness {
        factor_index: usize,
        ideal_norm: BigUint,
        generator_form: BinaryQuadraticForm,
    },
    /// A matched local oriented crater power could not be applied from the current node.
    LocalPower {
        factor_index: usize,
        source: OrientedCraterPowerActionError,
    },
}

impl fmt::Display for ClassGroupIsogenyActionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingLocalWitness {
                factor_index,
                ideal_norm,
                generator_form,
            } => {
                let human_factor_index = factor_index + 1;
                write!(
                    formatter,
                    "missing oriented local witness for factor {human_factor_index} with norm {ideal_norm} and form {generator_form:?}"
                )
            }
            Self::LocalPower {
                factor_index,
                source,
            } => write!(
                formatter,
                "could not apply oriented local power for factor {factor_index}: {source}"
            ),
        }
    }
}

impl std::error::Error for ClassGroupIsogenyActionError {}

/// One executed local factor in a staged class-group isogeny action.
///
/// The segment records a path obtained by applying one already-oriented local
/// crater power. It is geometric evidence inside a certified oriented crater,
/// not a proof that the supplied orientation is the arithmetic `𝔭` direction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassGroupIsogenyActionSegment {
    factor_index: usize,
    ideal: PrimeNormIdeal,
    generator_form: BinaryQuadraticForm,
    exponent: BigInt,
    start: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    path: Vec<IsogenyGraphNodeId>,
    direction_certification: CraterDirectionCertification,
}

impl ClassGroupIsogenyActionSegment {
    pub(crate) fn new(
        factor_index: usize,
        ideal: PrimeNormIdeal,
        generator_form: BinaryQuadraticForm,
        exponent: BigInt,
        path: Vec<IsogenyGraphNodeId>,
        direction_certification: CraterDirectionCertification,
    ) -> Self {
        let start = path
            .first()
            .copied()
            .expect("oriented local power reports have nonempty paths");
        let target = path
            .last()
            .copied()
            .expect("oriented local power reports have nonempty paths");

        Self {
            factor_index,
            ideal,
            generator_form,
            exponent,
            start,
            target,
            path,
            direction_certification,
        }
    }

    /// Returns the zero-based factor index in the source action plan.
    pub fn factor_index(&self) -> usize {
        self.factor_index
    }

    /// Returns the local prime-norm ideal attached to this segment.
    pub fn ideal(&self) -> &PrimeNormIdeal {
        &self.ideal
    }

    /// Returns the reduced form class attached to the local ideal.
    pub fn generator_form(&self) -> &BinaryQuadraticForm {
        &self.generator_form
    }

    /// Returns the exponent applied in the matched oriented crater.
    pub fn exponent(&self) -> &BigInt {
        &self.exponent
    }

    /// Returns the starting node for this local segment.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.start
    }

    /// Returns the final node reached by this local segment.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the local path, including both start and target.
    pub fn path(&self) -> &[IsogenyGraphNodeId] {
        &self.path
    }

    /// Returns how the local crater direction was certified.
    pub fn direction_certification(&self) -> CraterDirectionCertification {
        self.direction_certification
    }
}

/// Geometric execution report for an algebraic class-group action plan.
///
/// The report concatenates local oriented crater-power segments in the explicit
/// order chosen by the source plan. It remains a staged report: every segment
/// depends on a supplied local orientation witness, and no segment claims an
/// automatically inferred arithmetic orientation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClassGroupIsogenyActionReport {
    start: IsogenyGraphNodeId,
    target: IsogenyGraphNodeId,
    segments: Vec<ClassGroupIsogenyActionSegment>,
}

impl ClassGroupIsogenyActionReport {
    pub(crate) fn new(
        start: IsogenyGraphNodeId,
        segments: Vec<ClassGroupIsogenyActionSegment>,
    ) -> Self {
        let target = segments.last().map_or(start, |segment| segment.target());
        Self {
            start,
            target,
            segments,
        }
    }

    /// Returns the node where the staged action execution started.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.start
    }

    /// Returns the final node reached after all local segments.
    pub fn target(&self) -> IsogenyGraphNodeId {
        self.target
    }

    /// Returns the executed local segments in plan order.
    pub fn segments(&self) -> &[ClassGroupIsogenyActionSegment] {
        &self.segments
    }

    /// Returns whether the source plan had no nonzero local factors.
    pub fn is_trivial(&self) -> bool {
        self.segments.is_empty()
    }
}
