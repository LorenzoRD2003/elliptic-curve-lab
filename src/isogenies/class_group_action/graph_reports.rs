use std::{fmt, hash::Hash};

use crate::elliptic_curves::endomorphisms::quadratic_ideals::PrimeNormIdeal;
use crate::isogenies::{
    class_group_action::HorizontalIdealReport,
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
}
