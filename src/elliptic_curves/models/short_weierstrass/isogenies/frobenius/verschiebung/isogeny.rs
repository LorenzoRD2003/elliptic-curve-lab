use crate::elliptic_curves::short_weierstrass::{
    ShortWeierstrassCurve,
    function_fields::ShortWeierstrassFunction,
    isogenies::{
        frobenius::VerschiebungCertificate,
        frobenius::{AbsoluteFrobeniusIsogeny, FrobeniusLikeIsogeny},
        function_field_maps::ShortWeierstrassFunctionFieldMap,
    },
};
use crate::fields::traits::FiniteField;
use crate::isogenies::{
    error::{IsogenyError, VerschiebungError},
    traits::Isogeny,
};
use num_traits::ToPrimitive;

/// Function-field-side Verschiebung attached to one absolute Frobenius isogeny.
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
/// The current implementation stores only a pullback
/// `V* : F(E) -> F(E^(p))` and verifies these identities by composing
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

impl<F: FiniteField> VerschiebungIsogeny<F>
where
    F::Elem: PartialEq,
{
    /// Builds a Verschiebung pullback from an absolute Frobenius and a
    /// pullback map in the opposite direction.
    pub fn new(
        frobenius: AbsoluteFrobeniusIsogeny<F>,
        pullback: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<Self, IsogenyError> {
        if pullback.domain_curve() != frobenius.codomain()
            || pullback.codomain_curve() != frobenius.domain()
        {
            return Err(VerschiebungError::DomainCodomainMismatch.into());
        }

        Ok(Self {
            frobenius,
            pullback,
        })
    }

    /// Returns the absolute Frobenius this Verschiebung is meant to dualize.
    pub fn frobenius(&self) -> &AbsoluteFrobeniusIsogeny<F> {
        &self.frobenius
    }

    /// Returns the source curve `E^(p)` of `V : E^(p) -> E`.
    pub fn domain_curve(&self) -> &ShortWeierstrassCurve<F> {
        self.pullback.domain_curve()
    }

    /// Returns the target curve `E` of `V : E^(p) -> E`.
    pub fn codomain_curve(&self) -> &ShortWeierstrassCurve<F> {
        self.pullback.codomain_curve()
    }

    /// Returns the expected degree `p`.
    pub fn degree(&self) -> u128 {
        F::characteristic()
            .to_positive_biguint()
            .expect("finite fields have positive characteristic")
            .to_u128()
            .expect("Verschiebung degree should fit in u128 in the educational setting")
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
    ///
    /// Complexity: one pullback composition plus one map equality check.
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
            Err(VerschiebungError::LeftDualityViolation.into())
        }
    }

    /// Verifies `Frob_p ∘ V = [p]_{E^(p)}` by comparing pullbacks.
    ///
    /// Complexity: one pullback composition plus one map equality check.
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
            Err(VerschiebungError::RightDualityViolation.into())
        }
    }

    /// Verifies duality relations against supplied `[p]` pullback maps.
    ///
    /// Complexity: the sum of [`Self::verify_v_after_f_equals_p`] and
    /// [`Self::verify_f_after_v_equals_p`].
    pub fn verify_duality_relations(
        &self,
        multiplication_by_p_on_e: &ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: &ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<(), IsogenyError> {
        self.verify_v_after_f_equals_p(multiplication_by_p_on_e)?;
        self.verify_f_after_v_equals_p(multiplication_by_p_on_frobenius_twist)?;
        Ok(())
    }

    /// Packages this isogeny together with the two expected `[p]` pullbacks.
    ///
    /// Complexity: dominated by the immediate duality verification performed
    /// by [`VerschiebungCertificate::new`].
    pub fn certify(
        self,
        multiplication_by_p_on_e: ShortWeierstrassFunctionFieldMap<F>,
        multiplication_by_p_on_frobenius_twist: ShortWeierstrassFunctionFieldMap<F>,
    ) -> Result<VerschiebungCertificate<F>, IsogenyError> {
        VerschiebungCertificate::new(
            self,
            multiplication_by_p_on_e,
            multiplication_by_p_on_frobenius_twist,
        )
    }
}
