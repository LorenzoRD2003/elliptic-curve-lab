use crate::elliptic_curves::endomorphisms::binary_quadratic_forms::QuadraticClassGroup;
use crate::isogenies::class_group_action::{
    CraterDirectionCertification, CraterOrientationWitness, CraterOrientationWitnessError,
    LabeledCraterWalkReport, OrientedCraterClassOrderComparison,
    OrientedCraterClassOrderComparisonError, OrientedCraterPowerActionError,
    OrientedCraterPowerActionReport,
};
use crate::isogenies::graphs::IsogenyGraphNodeId;
use num_bigint::BigInt;

/// A labeled crater walk equipped with a user-supplied crater orientation.
///
/// This report records an explicit orientation witness for the certified crater
/// cycle. It still does not claim that the chosen positive direction was
/// inferred from arithmetic kernel data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrientedLabeledCraterWalkReport {
    labeled_walk: LabeledCraterWalkReport,
    orientation: CraterOrientationWitness,
}

impl OrientedLabeledCraterWalkReport {
    pub(crate) fn new(
        labeled_walk: LabeledCraterWalkReport,
        orientation: CraterOrientationWitness,
    ) -> Result<Self, CraterOrientationWitnessError> {
        let start = labeled_walk.walk().start();
        orientation
            .oriented_cycle_from(start)
            .ok_or(CraterOrientationWitnessError::MissingSuccessor { source: start })?;

        Ok(Self {
            labeled_walk,
            orientation,
        })
    }

    /// Returns the underlying graph-deterministic labeled walk report.
    pub fn labeled_walk(&self) -> &LabeledCraterWalkReport {
        &self.labeled_walk
    }

    /// Returns the user-supplied crater orientation witness.
    pub fn orientation(&self) -> &CraterOrientationWitness {
        &self.orientation
    }

    /// Returns how this oriented wrapper certifies the walk direction.
    pub fn direction_certification(&self) -> CraterDirectionCertification {
        CraterDirectionCertification::UserSuppliedArithmeticOrientation
    }

    /// Applies a local oriented crater power from `start`.
    ///
    /// The exponent is interpreted in the user-supplied crater orientation:
    /// positive exponents follow declared successors, negative exponents follow
    /// declared predecessors, and zero stays at `start`. This is a local
    /// oriented-crater operation, not a general class-group action.
    pub fn apply_power_from(
        &self,
        start: IsogenyGraphNodeId,
        exponent: BigInt,
    ) -> Result<OrientedCraterPowerActionReport, OrientedCraterPowerActionError> {
        OrientedCraterPowerActionReport::new(self.clone(), start, exponent)
    }

    /// Compares the generator class order with the oriented crater orbit length.
    ///
    /// This is a diagnostic for examples and staged class-group-action work.
    /// It computes the order of the reduced form class labeling the local ideal
    /// and compares it with the cycle length obtained by following the
    /// user-supplied crater orientation from `start`.
    ///
    /// Complexity: `O(h(D) · C + m)`, where `h(D)` is the class number of the
    /// quadratic order, `C` is the cost of one form composition/reduction, and
    /// `m` is the oriented crater cycle length.
    pub fn compare_generator_order(
        &self,
        class_group: &QuadraticClassGroup,
        start: IsogenyGraphNodeId,
    ) -> Result<OrientedCraterClassOrderComparison, OrientedCraterClassOrderComparisonError> {
        OrientedCraterClassOrderComparison::new(self.clone(), class_group, start)
    }
}
