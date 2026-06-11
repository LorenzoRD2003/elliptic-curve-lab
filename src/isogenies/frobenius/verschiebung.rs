use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassFunction};
use crate::fields::FiniteField;
use crate::isogenies::{
    AbsoluteFrobeniusIsogeny, FrobeniusLikeIsogeny, Isogeny, IsogenyError,
    ShortWeierstrassFunctionFieldMap, VerschiebungError,
};

/// Candidate Verschiebung attached to one absolute Frobenius isogeny.
///
/// This type is intentionally more modest than the explicit point-evaluable
/// isogenies in the crate.
///
/// Mathematically, for `Frob_p : E -> E^(p)`, a Verschiebung is a map
/// `V : E^(p) -> E` satisfying:
///
/// - `V ∘ Frob_p = [p]_E`,
/// - `Frob_p ∘ V = [p]_{E^(p)}`.
///
/// The current implementation stores only a candidate pullback
/// `V^* : F(E) -> F(E^(p))`  and verifies these identities by composing
/// function-field pullbacks.
///
/// It does **not** implement [`crate::isogenies::Isogeny`] yet, because the
/// repository does not currently expose a general procedure that turns a
/// pullback map on function fields into honest point evaluation `E^(p) -> E`.
#[derive(Clone, Debug)]
pub struct VerschiebungIsogeny<F: FiniteField> {
    frobenius: AbsoluteFrobeniusIsogeny<F>,
    pullback: ShortWeierstrassFunctionFieldMap<F>,
}

/// Certification object for a candidate Verschiebung. This packages:
///
/// - one absolute Frobenius `Frob_p : E -> E^(p)`
/// - one candidate pullback `V^* : F(E) -> F(E^(p))`
/// - one expected pullback for `[p]_E`
/// - one expected pullback for `[p]_{E^(p)}`
///
/// and turns the duality checks  `V ∘ Frob_p = [p]_E` and
/// `Frob_p ∘ V = [p]_{E^(p)}` into parameterless verification methods.
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

impl<F: FiniteField> VerschiebungIsogeny<F>
where
    F::Elem: PartialEq,
{
    /// Builds a candidate Verschiebung from an absolute Frobenius and a
    /// pullback map in the opposite direction.
    pub fn new(
        frobenius: AbsoluteFrobeniusIsogeny<F>,
        pullback: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<Self, IsogenyError> {
        if pullback.domain_curve() != frobenius.codomain()
            || pullback.codomain_curve() != frobenius.domain()
        {
            return Err(IsogenyError::Verschiebung(
                VerschiebungError::DomainCodomainMismatch,
            ));
        }

        Ok(Self {
            frobenius,
            pullback,
        })
    }

    /// Returns the absolute Frobenius this candidate is meant to dualize.
    pub fn frobenius(&self) -> &AbsoluteFrobeniusIsogeny<F> {
        &self.frobenius
    }

    /// Returns the source curve `E^(p)` of the candidate `V : E^(p) -> E`.
    pub fn domain_curve(&self) -> &ShortWeierstrassCurve<F> {
        self.pullback.domain_curve()
    }

    /// Returns the target curve `E` of the candidate `V : E^(p) -> E`.
    pub fn codomain_curve(&self) -> &ShortWeierstrassCurve<F> {
        self.pullback.codomain_curve()
    }

    /// Returns the expected degree `p`.
    pub fn degree(&self) -> u128 {
        u128::from(F::characteristic())
    }

    /// Returns the stored pullback `V^* : F(E) -> F(E^(p))`.
    pub fn as_function_field_map(&self) -> &ShortWeierstrassFunctionFieldMap<F> {
        &self.pullback
    }

    /// Returns the stored image `V^*(x)`.
    pub fn x_pullback(&self) -> &ShortWeierstrassFunction<F> {
        self.pullback.x_pullback()
    }

    /// Returns the stored image `V^*(y)`.
    pub fn y_pullback(&self) -> &ShortWeierstrassFunction<F> {
        self.pullback.y_pullback()
    }

    /// Verifies `V ∘ Frob_p = [p]_E` by comparing pullbacks.
    pub fn verify_v_after_f_equals_p(
        &self,
        multiplication_by_p_on_e: &ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<(), IsogenyError> {
        let composed = self
            .frobenius
            .as_function_field_map()
            .compose(&self.pullback)?;

        if &composed == multiplication_by_p_on_e {
            Ok(())
        } else {
            Err(IsogenyError::Verschiebung(
                VerschiebungError::LeftDualityViolation,
            ))
        }
    }

    /// Verifies `Frob_p ∘ V = [p]_{E^(p)}` by comparing pullbacks.
    pub fn verify_f_after_v_equals_p(
        &self,
        multiplication_by_p_on_frobenius_twist: &ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<(), IsogenyError> {
        let composed = self
            .pullback
            .compose(&self.frobenius.as_function_field_map())?;

        if &composed == multiplication_by_p_on_frobenius_twist {
            Ok(())
        } else {
            Err(IsogenyError::Verschiebung(
                VerschiebungError::RightDualityViolation,
            ))
        }
    }

    /// Verifies duality relations against supplied `[p]` pullback maps.
    pub fn verify_duality_relations(
        &self,
        multiplication_by_p_on_e: &ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: &ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<(), IsogenyError> {
        self.verify_v_after_f_equals_p(multiplication_by_p_on_e)?;
        self.verify_f_after_v_equals_p(multiplication_by_p_on_frobenius_twist)?;
        Ok(())
    }

    /// Packages this candidate together with the two expected `[p]` pullbacks.
    pub fn certify(
        self,
        multiplication_by_p_on_e: ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<VerschiebungCertificate<F>, IsogenyError> {
        let certificate = VerschiebungCertificate::new(
            self,
            multiplication_by_p_on_e,
            multiplication_by_p_on_frobenius_twist,
        )?;
        certificate.verify_duality_relations()?;
        Ok(certificate)
    }
}

impl<F: FiniteField> VerschiebungCertificate<F>
where
    F::Elem: PartialEq,
{
    /// Builds a certification object around a candidate Verschiebung.
    pub fn new(
        verschiebung: VerschiebungIsogeny<F>,
        multiplication_by_p_on_e: ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<Self, IsogenyError> {
        if multiplication_by_p_on_e.domain_curve() != verschiebung.frobenius.domain()
            || multiplication_by_p_on_e.codomain_curve() != verschiebung.frobenius.domain()
            || multiplication_by_p_on_frobenius_twist.domain_curve()
                != verschiebung.frobenius.codomain()
            || multiplication_by_p_on_frobenius_twist.codomain_curve()
                != verschiebung.frobenius.codomain()
        {
            return Err(IsogenyError::Verschiebung(
                VerschiebungError::DomainCodomainMismatch,
            ));
        }

        Ok(Self {
            verschiebung,
            multiplication_by_p_on_e,
            multiplication_by_p_on_frobenius_twist,
        })
    }

    /// Returns the certified candidate Verschiebung.
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
    pub fn verify_v_after_f_equals_p(&self) -> Result<(), IsogenyError> {
        self.verschiebung
            .verify_v_after_f_equals_p(&self.multiplication_by_p_on_e)
    }

    /// Verifies `Frob_p ∘ V = [p]_{E^(p)}`.
    pub fn verify_f_after_v_equals_p(&self) -> Result<(), IsogenyError> {
        self.verschiebung
            .verify_f_after_v_equals_p(&self.multiplication_by_p_on_frobenius_twist)
    }

    /// Verifies both duality relations.
    pub fn verify_duality_relations(&self) -> Result<(), IsogenyError> {
        self.verify_v_after_f_equals_p()?;
        self.verify_f_after_v_equals_p()?;
        Ok(())
    }
}
