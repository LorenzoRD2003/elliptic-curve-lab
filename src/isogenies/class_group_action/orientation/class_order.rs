use core::fmt;

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::{BinaryQuadraticForm, BinaryQuadraticFormError, QuadraticClassGroup},
    quadratic_ideals::PrimeNormIdeal,
};
use crate::isogenies::{
    class_group_action::OrientedLabeledCraterWalkReport, graphs::IsogenyGraphNodeId,
};

/// Outcome of comparing a class order with an oriented crater orbit length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrientedCraterClassOrderStatus {
    /// The algebraic class order equals the observed oriented orbit length.
    MatchesOrientedOrbit,
    /// The algebraic class order and observed oriented orbit length differ.
    OrientedOrbitLengthDiffers,
    /// The supplied orientation did not close an orbit from the requested start.
    OrbitDidNotClose,
}

/// Failure modes for comparing a generator class order with a crater orbit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OrientedCraterClassOrderComparisonError {
    /// The supplied class group has a different discriminant from the labeled ideal.
    ClassGroupDiscriminantMismatch,
    /// The generator form did not have a computable class order in that group.
    GeneratorOrder(BinaryQuadraticFormError),
    /// The requested start node is not part of the oriented crater cycle.
    StartOutsideOrientedCrater { start: IsogenyGraphNodeId },
}

impl fmt::Display for OrientedCraterClassOrderComparisonError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ClassGroupDiscriminantMismatch => write!(
                formatter,
                "class group discriminant does not match the labeled ideal order"
            ),
            Self::GeneratorOrder(error) => {
                write!(
                    formatter,
                    "could not compute generator class order: {error}"
                )
            }
            Self::StartOutsideOrientedCrater { start } => {
                write!(
                    formatter,
                    "start node {:?} is outside the oriented crater",
                    start
                )
            }
        }
    }
}

impl std::error::Error for OrientedCraterClassOrderComparisonError {}

/// Diagnostic comparison between an ideal/form class order and a crater orbit.
///
/// The report compares two already available quantities:
///
/// - the order of the reduced form class attached to the local ideal label;
/// - the length of the orbit obtained by following a user-supplied crater
///   orientation from a chosen start node.
///
/// Equality is useful evidence for examples, but it is not a proof that the
/// orientation is the arithmetic direction of `𝔭` rather than `\bar{𝔭}`, nor
/// that a general class-group action has been computed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrientedCraterClassOrderComparison {
    oriented_labeled_walk: OrientedLabeledCraterWalkReport,
    start: IsogenyGraphNodeId,
    class_order: usize,
    oriented_orbit_length: Option<usize>,
    status: OrientedCraterClassOrderStatus,
}

impl OrientedCraterClassOrderComparison {
    pub(super) fn new(
        oriented_labeled_walk: OrientedLabeledCraterWalkReport,
        class_group: &QuadraticClassGroup,
        start: IsogenyGraphNodeId,
    ) -> Result<Self, OrientedCraterClassOrderComparisonError> {
        let ideal_discriminant = oriented_labeled_walk
            .labeled_walk()
            .local_label()
            .ideal()
            .order()
            .discriminant();
        if class_group.discriminant() != ideal_discriminant {
            return Err(OrientedCraterClassOrderComparisonError::ClassGroupDiscriminantMismatch);
        }

        if !oriented_labeled_walk.orientation().contains_node(start) {
            return Err(
                OrientedCraterClassOrderComparisonError::StartOutsideOrientedCrater { start },
            );
        }

        let class_order = class_group
            .order_of_reduced_form(Self::generator_form_from(&oriented_labeled_walk))
            .map_err(OrientedCraterClassOrderComparisonError::GeneratorOrder)?;
        let oriented_orbit_length = oriented_labeled_walk
            .orientation()
            .oriented_cycle_from(start)
            .map(|cycle| cycle.len() - 1);
        let status = match oriented_orbit_length {
            Some(length) if length == class_order => {
                OrientedCraterClassOrderStatus::MatchesOrientedOrbit
            }
            Some(_) => OrientedCraterClassOrderStatus::OrientedOrbitLengthDiffers,
            None => OrientedCraterClassOrderStatus::OrbitDidNotClose,
        };

        Ok(Self {
            oriented_labeled_walk,
            start,
            class_order,
            oriented_orbit_length,
            status,
        })
    }

    /// Returns the oriented labeled walk being inspected.
    pub fn oriented_labeled_walk(&self) -> &OrientedLabeledCraterWalkReport {
        &self.oriented_labeled_walk
    }

    /// Returns the start node used for the oriented orbit.
    pub fn start(&self) -> IsogenyGraphNodeId {
        self.start
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
        Self::generator_form_from(&self.oriented_labeled_walk)
    }

    /// Returns the algebraic order of the generator form class.
    pub fn class_order(&self) -> usize {
        self.class_order
    }

    /// Returns the oriented crater orbit length, if the orbit closed.
    pub fn oriented_orbit_length(&self) -> Option<usize> {
        self.oriented_orbit_length
    }

    /// Returns the diagnostic comparison status.
    pub fn status(&self) -> OrientedCraterClassOrderStatus {
        self.status
    }

    fn generator_form_from(
        oriented_labeled_walk: &OrientedLabeledCraterWalkReport,
    ) -> &BinaryQuadraticForm {
        oriented_labeled_walk
            .labeled_walk()
            .form_label()
            .reduced_form()
    }
}
