use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::fields::FiniteField;
use crate::isogenies::{
    AbsoluteFrobeniusIsogeny, DifferentialPullbackReport, IsogenyError,
    ShortWeierstrassFunctionFieldMap, VerschiebungCertificate, VerschiebungIsogeny,
};
use crate::isogenies::{FrobeniusLikeIsogeny, Isogeny};

/// Report for the characteristic-`p` factorization `[p] = V ∘ Frob_p`.
///
/// The source of truth is intentionally small:
///
/// - the direct pullback `[p]^* : F(E) -> F(E)`
/// - a certificate containing the reconstructed Verschiebung and both duality
///   relations against `[p]_E` and `[p]_{E^(p)}`
///
/// Everything else in the report is derived from those two pieces.
#[derive(Debug)]
pub struct FrobeniusVerschiebungFactorizationReport<F: FiniteField> {
    multiplication_by_p_pullback: ShortWeierstrassFunctionFieldMap<F>,
    certificate: VerschiebungCertificate<F>,
}

impl<F: FiniteField> FrobeniusVerschiebungFactorizationReport<F>
where
    F::Elem: PartialEq,
{
    pub(crate) fn new(
        multiplication_by_p_pullback: ShortWeierstrassFunctionFieldMap<F>,
        certificate: VerschiebungCertificate<F>,
    ) -> Self {
        Self {
            multiplication_by_p_pullback,
            certificate,
        }
    }

    /// Returns the source curve `E`.
    pub fn curve(&self) -> &ShortWeierstrassCurve<F> {
        self.certificate.frobenius().domain()
    }

    /// Returns the characteristic `p`.
    pub fn scalar(&self) -> u64 {
        F::characteristic()
    }

    /// Returns the absolute Frobenius `Frob_p : E -> E^(p)`.
    pub fn frobenius(&self) -> &AbsoluteFrobeniusIsogeny<F> {
        self.certificate.frobenius()
    }

    /// Returns the direct pullback `[p]^* : F(E) -> F(E)`.
    pub fn multiplication_by_p_pullback(&self) -> &ShortWeierstrassFunctionFieldMap<F> {
        &self.multiplication_by_p_pullback
    }

    /// Returns the reconstructed Verschiebung `V : E^(p) -> E`.
    pub fn verschiebung(&self) -> &VerschiebungIsogeny<F> {
        self.certificate.verschiebung()
    }

    /// Returns the certified direct pullback `[p]^* : F(E^(p)) -> F(E^(p))`
    /// on the Frobenius twist.
    pub fn multiplication_by_p_on_frobenius_twist(&self) -> &ShortWeierstrassFunctionFieldMap<F> {
        self.certificate.multiplication_by_p_on_frobenius_twist()
    }

    /// Returns the stored certificate.
    pub fn certificate(&self) -> &VerschiebungCertificate<F> {
        &self.certificate
    }

    /// Recomputes the differential report for `Frob_p`.
    ///
    /// Complexity: the same as one call to
    /// [`ShortWeierstrassFunctionFieldMap::differential_pullback_report`] on
    /// the Frobenius pullback map.
    pub fn frobenius_differential_report(
        &self,
    ) -> Result<DifferentialPullbackReport<F>, IsogenyError> {
        self.frobenius().differential_pullback_report()
    }

    /// Recomputes the differential report for `V`.
    ///
    /// Complexity: the same as one call to
    /// [`ShortWeierstrassFunctionFieldMap::differential_pullback_report`] on
    /// the stored Verschiebung pullback map.
    pub fn verschiebung_differential_report(
        &self,
    ) -> Result<DifferentialPullbackReport<F>, IsogenyError> {
        self.verschiebung()
            .as_function_field_map()
            .differential_pullback_report()
    }

    /// Rechecks the two stored duality relations.
    ///
    /// Complexity: the same as
    /// [`VerschiebungCertificate::verify_duality_relations`], namely two
    /// pullback compositions plus their equality checks.
    pub fn verify(&self) -> Result<(), IsogenyError> {
        self.certificate.verify_duality_relations()
    }
}
