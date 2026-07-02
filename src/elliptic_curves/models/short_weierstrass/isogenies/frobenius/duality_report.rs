use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::isogenies::frobenius::VerschiebungCertificate,
};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    dual_report::{
        DegreeFactorizationSummary, DualIsogenyReport, DualIsogenySideSummary, DualityKind,
        KernelDescriptionSummary,
    },
    error::IsogenyError,
    traits::Isogeny,
};

impl<F: FiniteField> VerschiebungCertificate<F>
where
    F::Elem: PartialEq,
{
    /// Builds a structured report for the Frobenius/Verschiebung duality story.
    pub fn dual_report(
        &self,
    ) -> Result<DualIsogenyReport<ShortWeierstrassCurve<F>, ShortWeierstrassCurve<F>>, IsogenyError>
    {
        Ok(DualIsogenyReport::new(
            DualityKind::FrobeniusVerschiebung,
            DualIsogenySideSummary::new(
                self.frobenius().degree(),
                DegreeFactorizationSummary::from_degree_factorized_isogeny(self.frobenius()),
                KernelDescriptionSummary::from_kernel_description(&self.frobenius().kernel_description()),
            ),
            DualIsogenySideSummary::new(
                self.verschiebung().degree(),
                DegreeFactorizationSummary::unknown(),
                KernelDescriptionSummary {
                    total_degree: None,
                    reduced_degree: None,
                    infinitesimal_degree: None,
                    is_fully_reduced: false,
                    rational_point_count: None,
                    short_label: "kernel of Verschiebung is not yet modeled explicitly".to_string(),
                },
            ),
            self.verify_v_after_f_equals_p().is_ok(),
            self.verify_f_after_v_equals_p().is_ok(),
            vec![
                "This report describes the certified absolute-Frobenius/Verschiebung factorization.".to_string(),
                "The Frobenius side carries explicit degree factorization data.".to_string(),
                "The current crate does not yet model a fully explicit kernel description or separable/inseparable degree split for Verschiebung itself.".to_string(),
            ],
        ))
    }
}
