use super::VerschiebungIsogeny;
use crate::elliptic_curves::short_weierstrass::isogenies::{
    frobenius::AbsoluteFrobeniusIsogeny, function_field_maps::ShortWeierstrassFunctionFieldMap,
};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    error::{IsogenyError, VerschiebungError},
    traits::Isogeny,
};

/// Certification object for a validated Verschiebung. This packages:
///
/// - one absolute Frobenius `Frob_p : E -> E^(p)`
/// - one already shape-valid pullback `V* : F(E) -> F(E^(p))`
/// - one expected pullback for `[p]_E`
/// - one expected pullback for `[p]_{E^(p)}`
///
/// and turns the duality checks  `V ∘ Frob_p = [p]_E` and
/// `Frob_p ∘ V = [p]_{E^(p)}` into parameterless verification methods.
///
/// The constructor of this certificate verifies both relations immediately, so
/// every stored `VerschiebungIsogeny` is already safe in the sense that its
/// function-field data is compatible with the supplied direct `[p]` pullbacks.
///
/// This object is intentionally about certification data rather than point
/// evaluation. It still does **not** implement [`crate::isogenies::Isogeny`],
/// for the same reason as [`VerschiebungIsogeny`].
#[derive(Clone, Debug)]
pub struct VerschiebungCertificate<F: FiniteField> {
    verschiebung: VerschiebungIsogeny<F>,
    multiplication_by_p_on_e: ShortWeierstrassFunctionFieldMap<F>,
    multiplication_by_p_on_frobenius_twist: ShortWeierstrassFunctionFieldMap<F>,
}

impl<F: FiniteField> VerschiebungCertificate<F>
where
    F::Elem: PartialEq,
{
    /// Builds a certification object around a validated Verschiebung pullback.
    ///
    /// Complexity: `Θ(1)` curve-compatibility checks plus two pullback
    /// compositions and equality checks coming from the immediate duality
    /// verification.
    pub fn new(
        verschiebung: VerschiebungIsogeny<F>,
        multiplication_by_p_on_e: ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<Self, IsogenyError> {
        if multiplication_by_p_on_e.domain_curve() != verschiebung.frobenius().domain()
            || multiplication_by_p_on_e.codomain_curve() != verschiebung.frobenius().domain()
            || multiplication_by_p_on_frobenius_twist.domain_curve()
                != verschiebung.frobenius().codomain()
            || multiplication_by_p_on_frobenius_twist.codomain_curve()
                != verschiebung.frobenius().codomain()
        {
            return Err(VerschiebungError::DomainCodomainMismatch.into());
        }

        let certificate = Self {
            verschiebung,
            multiplication_by_p_on_e,
            multiplication_by_p_on_frobenius_twist,
        };
        certificate.verify_duality_relations()?;
        Ok(certificate)
    }

    /// Returns the certified Verschiebung isogeny.
    pub fn verschiebung(&self) -> &VerschiebungIsogeny<F> {
        &self.verschiebung
    }

    /// Returns the stored absolute Frobenius.
    pub fn frobenius(&self) -> &AbsoluteFrobeniusIsogeny<F> {
        self.verschiebung.frobenius()
    }

    /// Returns the stored expected pullback for `[p]_E`.
    pub fn multiplication_by_p_on_e(&self) -> &ShortWeierstrassFunctionFieldMap<F> {
        &self.multiplication_by_p_on_e
    }

    /// Returns the stored expected pullback for `[p]_{E^(p)}`.
    pub fn multiplication_by_p_on_frobenius_twist(&self) -> &ShortWeierstrassFunctionFieldMap<F> {
        &self.multiplication_by_p_on_frobenius_twist
    }

    /// Verifies `V ∘ Frob_p = [p]_E`.
    ///
    /// Complexity: one pullback composition plus one map equality check.
    pub fn verify_v_after_f_equals_p(&self) -> Result<(), IsogenyError> {
        self.verschiebung
            .verify_v_after_f_equals_p(&self.multiplication_by_p_on_e)
    }

    /// Verifies `Frob_p ∘ V = [p]_{E^(p)}`.
    ///
    /// Complexity: one pullback composition plus one map equality check.
    pub fn verify_f_after_v_equals_p(&self) -> Result<(), IsogenyError> {
        self.verschiebung
            .verify_f_after_v_equals_p(&self.multiplication_by_p_on_frobenius_twist)
    }

    /// Verifies both duality relations.
    ///
    /// Complexity: the sum of [`Self::verify_v_after_f_equals_p`] and
    /// [`Self::verify_f_after_v_equals_p`].
    pub fn verify_duality_relations(&self) -> Result<(), IsogenyError> {
        self.verify_v_after_f_equals_p()?;
        self.verify_f_after_v_equals_p()?;
        Ok(())
    }
}
