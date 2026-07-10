use std::{fmt, hash::Hash};

use crate::elliptic_curves::endomorphisms::{
    binary_quadratic_forms::QuadraticClassGroup, quadratic_ideals::PrimeNormIdeal,
};
use crate::isogenies::graphs::IsogenyGraphNodeId;
use crate::isogenies::{
    class_group_action::{
        CraterWalkReport, HorizontalIdealReport, LabeledCraterWalkError, LabeledCraterWalkReport,
    },
    graphs::{GraphCurveModel, IsogenyGraph, endomorphisms::VolcanoSearchError},
};

impl<C: GraphCurveModel> IsogenyGraph<C>
where
    C::Point: Clone + Eq + Hash,
    C::IsomorphismWitness: Clone + fmt::Debug,
{
    /// Builds horizontal ideal reports for the supplied prime-norm ideal.
    ///
    /// The ideal norm selects the local volcano prime `ℓ`; the graph then
    /// computes its crater report at `ℓ` and annotates the reported horizontal
    /// edges with the same ideal.
    ///
    /// Complexity: one crater-report construction for `ℓ`, followed by linear
    /// report assembly in the number of crater-horizontal edge reports.
    pub fn horizontal_ideal_reports(
        &self,
        ideal: PrimeNormIdeal,
    ) -> Result<Vec<HorizontalIdealReport>, VolcanoSearchError> {
        let crater = self.volcano_crater_report(ideal.norm())?;
        Ok(HorizontalIdealReport::for_crater_report(&crater, ideal))
    }

    /// Builds a deterministic crater walk labeled by the supplied ideal.
    ///
    /// The ideal norm selects the local volcano prime `ℓ`. The graph computes
    /// its crater report at `ℓ`, starts at `start`, and follows certified
    /// horizontal crater edges in a deterministic graph order.
    ///
    /// The returned report records the visited path and the cycle length when
    /// the walk closes back at `start`. See [`CraterWalkReport`] for the local
    /// direction convention.
    ///
    /// Complexity: one crater-report construction for `ℓ`, plus a linear scan
    /// of crater-horizontal edge reports.
    pub fn crater_walk_report(
        &self,
        ideal: PrimeNormIdeal,
        start: IsogenyGraphNodeId,
    ) -> Result<CraterWalkReport, VolcanoSearchError> {
        let crater = self.volcano_crater_report(ideal.norm())?;

        Ok(CraterWalkReport::from_crater_report(&crater, ideal, start))
    }

    /// Builds a deterministic crater walk with its ideal/form-class labels.
    ///
    /// This is a compatibility-and-visualization report, not a certified
    /// class-group action. The ideal norm selects the local volcano prime `ℓ`;
    /// the class group supplies the discriminant used to validate the local
    /// label; and the walk direction remains graph-deterministic.
    ///
    /// Complexity: one crater-report construction for `ℓ`, one local
    /// ideal/class-group compatibility check, one ideal-to-form conversion, and
    /// one deterministic crater walk through certified crater edges.
    pub fn labeled_crater_walk_report(
        &self,
        class_group: &QuadraticClassGroup,
        ideal: PrimeNormIdeal,
        start: IsogenyGraphNodeId,
    ) -> Result<LabeledCraterWalkReport, LabeledCraterWalkError> {
        let crater = self.volcano_crater_report(ideal.norm())?;

        Ok(LabeledCraterWalkReport::from_crater_report(
            &crater,
            class_group,
            ideal,
            start,
        )?)
    }
}
