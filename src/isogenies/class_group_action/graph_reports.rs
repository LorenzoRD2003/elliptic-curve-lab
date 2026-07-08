use std::{fmt, hash::Hash};

use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;
use crate::isogenies::graphs::IsogenyGraphNodeId;
use crate::isogenies::{
    class_group_action::{CraterWalkReport, HorizontalIdealReport},
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
}
