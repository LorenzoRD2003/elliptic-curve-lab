use crate::fields::traits::*;
use std::hash::Hash;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::isogenies::{DualVeluIsogeny, VeluIsogeny},
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    dual_report::{
        DegreeFactorizationSummary, DualIsogenyReport, DualIsogenySideSummary, DualityKind,
        KernelDescriptionSummary,
    },
    error::IsogenyError,
    traits::Isogeny,
};

impl<F: Field + EnumerableFiniteField + SqrtField + Clone> DualVeluIsogeny<F>
where
    F::Elem: Clone + Eq + Hash,
{
    /// Builds a structured report for one classical separable dual pair.
    pub fn dual_report(
        phi: &VeluIsogeny<ShortWeierstrassCurve<F>>,
        dual: &Self,
    ) -> Result<DualIsogenyReport<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>, IsogenyError>
    {
        Ok(DualIsogenyReport::new(
            DualityKind::SeparableClassical,
            DualIsogenySideSummary::new(
                phi.degree(),
                DegreeFactorizationSummary::from_degree_factorized_isogeny(phi),
                KernelDescriptionSummary::from_kernel_description(&phi.kernel_description()),
            ),
            DualIsogenySideSummary::new(
                dual.degree(),
                DegreeFactorizationSummary::from_degree_factorized_isogeny(dual),
                KernelDescriptionSummary::from_kernel_description(&dual.kernel_description()),
            ),
            dual.verify_left_dual_relation(phi).is_ok(),
            dual.verify_right_dual_relation(phi).is_ok(),
            vec![
                "This report describes the classical small-field duality search for a separable Velu isogeny.".to_string(),
                "In the current educational Velu layer, the dual pair is modeled as separable with inseparable degree 1.".to_string(),
            ],
        ))
    }
}
