use core::fmt;

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::ToPrimitive;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::BinaryQuadraticForm, quadratic_ideals::PrimeNormIdeal,
};
use crate::isogenies::{
    class_group_action::OrientedLabeledCraterWalkReport, graphs::IsogenyGraphNodeId,
};

/// Failure modes for applying a local oriented crater power.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OrientedCraterPowerActionError {
    /// The requested start node is not part of the oriented crater cycle.
    StartOutsideOrientedCrater { start: IsogenyGraphNodeId },
    /// The orientation witness did not provide a required successor/predecessor.
    MissingSuccessor { source: IsogenyGraphNodeId },
    /// The computed path was unexpectedly empty.
    EmptyPath,
}

impl fmt::Display for OrientedCraterPowerActionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StartOutsideOrientedCrater { start } => {
                write!(
                    formatter,
                    "start node {:?} is outside the oriented crater",
                    start
                )
            }
            Self::MissingSuccessor { source } => {
                write!(formatter, "orientation has no next step from {:?}", source)
            }
            Self::EmptyPath => write!(formatter, "oriented crater power produced an empty path"),
        }
    }
}

impl std::error::Error for OrientedCraterPowerActionError {}

/// Local power of a prime-ideal label inside a user-oriented crater.
///
/// The report records a path obtained by walking an exponent number of steps in
/// the user-supplied crater orientation. It models powers of the selected local
/// ideal only under that explicit convention; it is not a general class-group
/// action constructor.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrientedCraterPowerActionReport {
    oriented_labeled_walk: OrientedLabeledCraterWalkReport,
    exponent: BigInt,
    path: Vec<IsogenyGraphNodeId>,
}

impl OrientedCraterPowerActionReport {
    pub(super) fn new(
        oriented_labeled_walk: OrientedLabeledCraterWalkReport,
        start: IsogenyGraphNodeId,
        exponent: BigInt,
    ) -> Result<Self, OrientedCraterPowerActionError> {
        if !oriented_labeled_walk.orientation().contains_node(start) {
            return Err(OrientedCraterPowerActionError::StartOutsideOrientedCrater { start });
        }

        let cycle_length = oriented_labeled_walk
            .orientation()
            .oriented_cycle_from(start)
            .ok_or(OrientedCraterPowerActionError::MissingSuccessor { source: start })?
            .len()
            - 1;
        let step_count = reduced_step_count(&exponent, cycle_length);
        let mut path = vec![start];
        let mut current = start;

        for _ in 0..step_count {
            current = match exponent.sign() {
                Sign::Minus => oriented_labeled_walk
                    .orientation()
                    .predecessor(current)
                    .ok_or(OrientedCraterPowerActionError::MissingSuccessor { source: current })?,
                Sign::NoSign | Sign::Plus => oriented_labeled_walk
                    .orientation()
                    .successor(current)
                    .ok_or(OrientedCraterPowerActionError::MissingSuccessor { source: current })?,
            };
            path.push(current);
        }

        if path.is_empty() {
            return Err(OrientedCraterPowerActionError::EmptyPath);
        }

        Ok(Self {
            oriented_labeled_walk,
            exponent,
            path,
        })
    }

    /// Returns the oriented labeled walk that supplies the local generator.
    pub fn oriented_labeled_walk(&self) -> &OrientedLabeledCraterWalkReport {
        &self.oriented_labeled_walk
    }

    /// Returns the exponent applied in the user-supplied orientation.
    pub fn exponent(&self) -> &BigInt {
        &self.exponent
    }

    /// Returns the recorded oriented path.
    pub fn path(&self) -> &[IsogenyGraphNodeId] {
        &self.path
    }

    /// Returns the final node reached by this local oriented power.
    pub fn target(&self) -> IsogenyGraphNodeId {
        *self
            .path
            .last()
            .expect("oriented crater power paths are nonempty")
    }

    /// Returns the local prime-norm ideal used as generator label.
    pub fn generator_ideal(&self) -> &PrimeNormIdeal {
        self.oriented_labeled_walk
            .labeled_walk()
            .local_label()
            .ideal()
    }

    /// Returns the reduced form class labeling the local generator.
    pub fn generator_form(&self) -> &BinaryQuadraticForm {
        self.oriented_labeled_walk
            .labeled_walk()
            .form_label()
            .reduced_form()
    }
}

fn reduced_step_count(exponent: &BigInt, cycle_length: usize) -> usize {
    if cycle_length == 0 {
        return 0;
    }
    let modulus = BigUint::from(cycle_length);
    (exponent.magnitude() % modulus)
        .to_usize()
        .expect("remainder modulo a usize should fit in usize")
}
