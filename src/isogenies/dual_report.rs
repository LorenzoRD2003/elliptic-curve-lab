use core::marker::PhantomData;
use std::hash::Hash;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::CurveModel;
use crate::fields::{EnumerableFiniteField, Field, FiniteField, SqrtField};
use crate::isogenies::{
    DegreeFactorizedIsogeny, DualVeluIsogeny, Isogeny, IsogenyError, KernelDescription,
    VeluIsogeny, VerschiebungCertificate, verify_left_dual_relation, verify_right_dual_relation,
};

/// High-level classification of the duality story represented by one report.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DualityKind {
    /// Classical duality for a separable small-field Vélu isogeny.
    SeparableClassical,
    /// Frobenius/Verschiebung duality in characteristic `p`.
    FrobeniusVerschiebung,
    /// Placeholder bucket for future partially modeled or mixed duality stories.
    MixedOrPartial,
}

/// Compact summary of separable/inseparable degree data.
///
/// When both entries are present, the intended interpretation is
///
/// `total degree = separable_degree * inseparable_degree`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DegreeFactorizationSummary {
    separable_degree: Option<u128>,
    inseparable_degree: Option<u128>,
}

/// Compact summary of the currently available kernel description.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KernelDescriptionSummary {
    total_degree: Option<usize>,
    reduced_degree: Option<usize>,
    infinitesimal_degree: Option<usize>,
    is_fully_reduced: bool,
    rational_point_count: Option<usize>,
    short_label: String,
}

/// One side of a duality report, either `phi` or `phi_hat`.
///
/// This groups together the ordinary degree, any finer
/// separable/inseparable factorization data, and the currently available
/// kernel summary for that same isogeny. Grouping them avoids parallel
/// `phi_*` / `dual_*` fields inside [`DualIsogenyReport`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DualIsogenySideSummary {
    degree: usize,
    degree_factorization: DegreeFactorizationSummary,
    kernel_summary: KernelDescriptionSummary,
}

/// Structured duality report for two related isogeny objects.
///
/// This report is intentionally lightweight. It does not duplicate full curve
/// or map data. Instead, it records:
///
/// - the kind of duality story currently being modeled
/// - one summary for `phi`
/// - one summary for `phi_hat`
/// - whether the left and right duality relations currently verify
/// - free-form notes explaining the current mathematical scope
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DualIsogenyReport<Domain: CurveModel, Codomain: CurveModel> {
    duality_kind: DualityKind,
    phi: DualIsogenySideSummary,
    dual: DualIsogenySideSummary,
    left_relation_holds: bool,
    right_relation_holds: bool,
    notes: Vec<String>,
    marker: PhantomData<(Domain, Codomain)>,
}

impl DegreeFactorizationSummary {
    /// Builds a degree summary with no finer separable/inseparable data.
    pub fn unknown() -> Self {
        Self {
            separable_degree: None,
            inseparable_degree: None,
        }
    }

    /// Builds a degree summary from explicit optional factors.
    pub fn new(separable_degree: Option<u128>, inseparable_degree: Option<u128>) -> Self {
        Self {
            separable_degree,
            inseparable_degree,
        }
    }

    /// Builds a degree summary from an isogeny that already models the
    /// separable/inseparable factorization explicitly.
    pub fn from_degree_factorized_isogeny<I, D: CurveModel, C: CurveModel>(isogeny: &I) -> Self
    where
        I: DegreeFactorizedIsogeny<D, C>,
    {
        Self::new(
            Some(isogeny.separable_degree()),
            Some(isogeny.inseparable_degree()),
        )
    }

    /// Returns the separable-degree factor when known.
    pub fn separable_degree(&self) -> Option<u128> {
        self.separable_degree
    }

    /// Returns the inseparable-degree factor when known.
    pub fn inseparable_degree(&self) -> Option<u128> {
        self.inseparable_degree
    }
}

impl KernelDescriptionSummary {
    /// Builds a compact summary from a full kernel description.
    pub fn from_kernel_description<C: CurveModel>(description: &KernelDescription<C>) -> Self {
        let short_label = match description {
            KernelDescription::Reduced(_) => "reduced kernel".to_string(),
            KernelDescription::NonReduced(nonreduced) => {
                format!("nonreduced kernel ({})", nonreduced.label())
            }
            KernelDescription::Mixed(mixed) => match mixed.label() {
                Some(label) => format!("mixed kernel ({label})"),
                None => "mixed kernel".to_string(),
            },
            KernelDescription::Unknown => "unknown kernel description".to_string(),
        };

        Self {
            total_degree: description.degree(),
            reduced_degree: description.reduced_degree(),
            infinitesimal_degree: description.infinitesimal_degree(),
            is_fully_reduced: description.is_fully_reduced(),
            rational_point_count: description.rational_points().map(<[C::Point]>::len),
            short_label,
        }
    }

    /// Returns the currently known total kernel degree.
    pub fn total_degree(&self) -> Option<usize> {
        self.total_degree
    }

    /// Returns the currently known reduced kernel degree.
    pub fn reduced_degree(&self) -> Option<usize> {
        self.reduced_degree
    }

    /// Returns the currently known infinitesimal kernel degree.
    pub fn infinitesimal_degree(&self) -> Option<usize> {
        self.infinitesimal_degree
    }

    /// Returns whether the current kernel description is fully reduced.
    pub fn is_fully_reduced(&self) -> bool {
        self.is_fully_reduced
    }

    /// Returns how many explicit rational points are currently visible in the
    /// kernel description, when that notion makes sense for the current case.
    pub fn rational_point_count(&self) -> Option<usize> {
        self.rational_point_count
    }

    /// Returns the short human-readable kernel label used by reports and
    /// visualization helpers.
    pub fn short_label(&self) -> &str {
        &self.short_label
    }
}

impl DualIsogenySideSummary {
    fn new(
        degree: usize,
        degree_factorization: DegreeFactorizationSummary,
        kernel_summary: KernelDescriptionSummary,
    ) -> Self {
        Self {
            degree,
            degree_factorization,
            kernel_summary,
        }
    }

    /// Returns the total degree recorded for this side.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Returns the separable/inseparable degree summary for this side.
    pub fn degree_factorization(&self) -> &DegreeFactorizationSummary {
        &self.degree_factorization
    }

    /// Returns the current kernel summary for this side.
    pub fn kernel_summary(&self) -> &KernelDescriptionSummary {
        &self.kernel_summary
    }
}

impl<Domain: CurveModel, Codomain: CurveModel> DualIsogenyReport<Domain, Codomain> {
    fn new(
        duality_kind: DualityKind,
        phi: DualIsogenySideSummary,
        dual: DualIsogenySideSummary,
        left_relation_holds: bool,
        right_relation_holds: bool,
        notes: Vec<String>,
    ) -> Self {
        Self {
            duality_kind,
            phi,
            dual,
            left_relation_holds,
            right_relation_holds,
            notes,
            marker: PhantomData,
        }
    }

    /// Returns the high-level duality kind modeled by this report.
    pub fn duality_kind(&self) -> DualityKind {
        self.duality_kind
    }

    /// Returns the grouped summary for `phi`.
    pub fn phi(&self) -> &DualIsogenySideSummary {
        &self.phi
    }

    /// Returns the grouped summary for `phi_hat`.
    pub fn dual(&self) -> &DualIsogenySideSummary {
        &self.dual
    }

    /// Returns `deg(phi)`.
    pub fn phi_degree(&self) -> usize {
        self.phi.degree()
    }

    /// Returns `deg(phi_hat)`.
    pub fn dual_degree(&self) -> usize {
        self.dual.degree()
    }

    /// Returns the separable/inseparable degree summary for `phi`.
    pub fn phi_degree_factorization(&self) -> &DegreeFactorizationSummary {
        self.phi.degree_factorization()
    }

    /// Returns the separable/inseparable degree summary for `phi_hat`.
    pub fn dual_degree_factorization(&self) -> &DegreeFactorizationSummary {
        self.dual.degree_factorization()
    }

    /// Returns the kernel summary for `phi`.
    pub fn phi_kernel_summary(&self) -> &KernelDescriptionSummary {
        self.phi.kernel_summary()
    }

    /// Returns the kernel summary for `phi_hat`.
    pub fn dual_kernel_summary(&self) -> &KernelDescriptionSummary {
        self.dual.kernel_summary()
    }

    /// Returns whether the left duality relation currently verifies.
    pub fn left_relation_holds(&self) -> bool {
        self.left_relation_holds
    }

    /// Returns whether the right duality relation currently verifies.
    pub fn right_relation_holds(&self) -> bool {
        self.right_relation_holds
    }

    /// Returns explanatory notes about the current mathematical scope.
    pub fn notes(&self) -> &[String] {
        &self.notes
    }
}

impl<F> DualVeluIsogeny<F>
where
    F: Field + EnumerableFiniteField + SqrtField + Clone,
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
            verify_left_dual_relation(phi, dual).is_ok(),
            verify_right_dual_relation(phi, dual).is_ok(),
            vec![
                "This report describes the classical small-field duality search for a separable Velu isogeny.".to_string(),
                "In the current educational Velu layer, the dual pair is modeled as separable with inseparable degree 1.".to_string(),
            ],
        ))
    }
}

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
                usize::try_from(self.verschiebung().degree())
                    .expect("educational Verschiebung degree should fit into usize"),
                DegreeFactorizationSummary::unknown(),
                KernelDescriptionSummary {
                    total_degree: usize::try_from(self.verschiebung().degree()).ok(),
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

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{AffineCurveModel, ShortWeierstrassCurve};
    use crate::fields::{Field, Fp};
    use crate::isogenies::{
        AbsoluteFrobeniusIsogeny, DualIsogenyReport, DualVeluIsogeny, DualityKind,
        FrobeniusLikeIsogeny, Isogeny, VeluIsogeny, VerschiebungCertificate, VerschiebungIsogeny,
    };

    type F29 = Fp<29>;
    type F41 = Fp<41>;
    type Curve29 = ShortWeierstrassCurve<F29>;
    type Curve41 = ShortWeierstrassCurve<F41>;

    fn curve_f29() -> Curve29 {
        Curve29::new(F29::from_i64(2), F29::from_i64(2)).expect("valid curve")
    }

    fn curve_f41() -> Curve41 {
        Curve41::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    fn degree_three_phi() -> VeluIsogeny<Curve29> {
        let curve = curve_f29();
        let generator = curve
            .point(F29::from_i64(10), F29::from_i64(23))
            .expect("sample generator should lie on the curve");
        VeluIsogeny::from_generator(curve, generator).expect("sample Velu isogeny should build")
    }

    #[test]
    fn classical_dual_report_records_basic_duality_data() {
        let phi = degree_three_phi();
        let dual = phi.find_dual_exhaustively().expect("dual should be found");
        let report: DualIsogenyReport<Curve29, Curve29> =
            DualVeluIsogeny::dual_report(&phi, &dual).expect("report should build");

        assert_eq!(report.duality_kind(), DualityKind::SeparableClassical);
        assert_eq!(report.phi_degree(), 3);
        assert_eq!(report.dual_degree(), 3);
        assert!(report.left_relation_holds());
        assert!(report.right_relation_holds());
        assert_eq!(report.phi_kernel_summary().total_degree(), Some(3));
        assert_eq!(report.dual_kernel_summary().total_degree(), Some(3));
        assert_eq!(
            report.phi_degree_factorization().separable_degree(),
            Some(3)
        );
        assert_eq!(
            report.phi_degree_factorization().inseparable_degree(),
            Some(1)
        );
        assert_eq!(
            report.dual_degree_factorization().separable_degree(),
            Some(3)
        );
        assert_eq!(
            report.dual_degree_factorization().inseparable_degree(),
            Some(1)
        );
    }

    #[test]
    fn frobenius_verschiebung_report_records_known_frobenius_data_and_unknown_v_data() {
        let curve = curve_f41();
        let frobenius =
            AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
        let candidate_v = crate::isogenies::ShortWeierstrassFunctionFieldMap::new(
            frobenius.codomain().clone(),
            frobenius.domain().clone(),
            frobenius
                .as_function_field_map()
                .codomain_function_field()
                .x(),
            frobenius
                .as_function_field_map()
                .codomain_function_field()
                .y(),
        )
        .expect("identity candidate on the twist should validate");
        let verschiebung = VerschiebungIsogeny::new(frobenius.clone(), candidate_v)
            .expect("candidate should build");
        let expected_left = frobenius
            .as_function_field_map()
            .compose(verschiebung.as_function_field_map())
            .expect("left composition should build");
        let expected_right = verschiebung
            .as_function_field_map()
            .compose(&frobenius.as_function_field_map())
            .expect("right composition should build");
        let certificate = VerschiebungCertificate::new(verschiebung, expected_left, expected_right)
            .expect("certificate should build");
        let report = certificate.dual_report().expect("report should build");

        assert_eq!(report.duality_kind(), DualityKind::FrobeniusVerschiebung);
        assert_eq!(report.phi_degree(), 41);
        assert_eq!(report.dual_degree(), 41);
        assert_eq!(
            report.phi_degree_factorization().separable_degree(),
            Some(1)
        );
        assert_eq!(
            report.phi_degree_factorization().inseparable_degree(),
            Some(41)
        );
        assert!(
            report
                .dual_degree_factorization()
                .separable_degree()
                .is_none()
        );
        assert!(report.left_relation_holds());
        assert!(report.right_relation_holds());
        assert!(
            report
                .dual_kernel_summary()
                .short_label()
                .contains("not yet modeled")
        );
    }
}
