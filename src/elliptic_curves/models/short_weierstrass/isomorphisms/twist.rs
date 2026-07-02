use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::quadratic_twist::QuadraticTwistFrobeniusRelation,
    short_weierstrass::isomorphisms::{CurveIsomorphismError, ShortWeierstrassIsomorphism},
    traits::{CurveIsomorphism, FrobeniusTraceCurveModel},
};
use crate::fields::traits::*;
use crate::fields::{
    FieldError,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    polynomial_field::PolynomialModulus,
    traits::{EnumerableFiniteField, FiniteField, SqrtField},
};

/// Whether a quadratic twist is trivial or genuinely quadratic over the
/// current base field.
///
/// The generic short-Weierstrass case `j != 0, 1728` admits only the
/// geometric automorphisms `{±1}`, so a twist by `d` is already base-field
/// trivial exactly when the usual square-root witness exists.
///
/// The exceptional families `j = 1728` (`b = 0`) and `j = 0` (`a = 0`) admit
/// extra geometric automorphisms. In particular, for `j = 1728` a non-square
/// twist factor can still produce a base-field-trivial twist. The current
/// implementation certifies that extra `j = 1728` path explicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TwistKind {
    Trivial,
    Quadratic,
}

/// Quadratic-twist package for a short-Weierstrass curve.
///
/// The primary parameter is the twist factor `d`. If `d` is a square in the base field,
/// then the twist is already trivial over that field and a base-field isomorphism
/// exists. If `d` is not a square, then the stored twist is typically only isomorphic
/// over a field extension containing `sqrt(d)`.
pub struct ShortWeierstrassQuadraticTwist<F: Field> {
    original: ShortWeierstrassCurve<F>,
    twist: ShortWeierstrassCurve<F>,
    d: F::Elem,
}

/// Backward-compatible alias for the current short-Weierstrass quadratic twist
/// package.
pub type ShortWeierstrassTwist<F> = ShortWeierstrassQuadraticTwist<F>;

impl<F: Field> ShortWeierstrassQuadraticTwist<F> {
    /// Builds the quadratic twist package determined by the factor `d`.
    pub fn new(
        original: ShortWeierstrassCurve<F>,
        d: F::Elem,
    ) -> Result<Self, CurveIsomorphismError> {
        let twist = original.quadratic_twist(d.clone())?;
        Ok(Self { original, twist, d })
    }

    /// Returns the original curve `E`.
    pub fn original(&self) -> &ShortWeierstrassCurve<F> {
        &self.original
    }

    /// Returns the twisted curve `E^(d)`.
    pub fn twist(&self) -> &ShortWeierstrassCurve<F> {
        &self.twist
    }

    /// Returns the twist factor `d`.
    pub fn factor(&self) -> &F::Elem {
        &self.d
    }

    fn original_curve_copy(&self) -> ShortWeierstrassCurve<F> {
        ShortWeierstrassCurve::new(self.original.a().clone(), self.original.b().clone())
            .expect("stored original curve should stay non-singular")
    }

    fn expected_quadratic_extension_modulus(
        &self,
    ) -> Result<PolynomialModulus<F>, CurveIsomorphismError> {
        PolynomialModulus::new(vec![F::neg(&self.d), F::zero(), F::one()]).map_err(Into::into)
    }

    fn lift_curve_to_extension<S: ExtensionFieldSpec<Base = F>>(
        curve: &ShortWeierstrassCurve<F>,
    ) -> Result<ShortWeierstrassCurve<ExtensionField<S>>, CurveIsomorphismError> {
        ShortWeierstrassCurve::<ExtensionField<S>>::new(
            ExtensionField::from_base(curve.a().clone()),
            ExtensionField::from_base(curve.b().clone()),
        )
        .map_err(Into::into)
    }
}

impl<F: SqrtField> ShortWeierstrassQuadraticTwist<F> {
    /// Returns whether the stored twist is trivial or genuinely quadratic over
    /// the current base field.
    pub fn kind(&self) -> TwistKind {
        if self.certified_base_field_scaling_factor().is_some() {
            TwistKind::Trivial
        } else {
            TwistKind::Quadratic
        }
    }

    /// Returns one certified base-field isomorphism from `E` to `E^(d)` when
    /// the current field backend can witness it directly.
    pub fn base_field_isomorphism(&self) -> Option<ShortWeierstrassIsomorphism<F>> {
        let scaling_factor = self.certified_base_field_scaling_factor()?;
        ShortWeierstrassIsomorphism::new(self.original_curve_copy(), scaling_factor).ok()
    }

    /// Returns the canonical scaling isomorphism over a genuine quadratic
    /// extension presented as `F[x] / (x^2 - d)`.
    pub fn isomorphism_over_quadratic_extension<S: ExtensionFieldSpec<Base = F>>(
        &self,
    ) -> Result<ShortWeierstrassIsomorphism<ExtensionField<S>>, CurveIsomorphismError> {
        if self.kind() == TwistKind::Trivial {
            return Err(FieldError::NonIrreduciblePolynomial.into());
        }

        let expected_modulus = self.expected_quadratic_extension_modulus()?;
        if S::defining_modulus() != expected_modulus {
            return Err(FieldError::IncompatibleFieldParameters.into());
        }

        ExtensionField::<S>::check_structure()?;

        let lifted_domain = Self::lift_curve_to_extension(&self.original)?;
        let lifted_twist = Self::lift_curve_to_extension(&self.twist)?;
        let u = ExtensionField::<S>::element(vec![F::zero(), F::one()]);
        let isomorphism = ShortWeierstrassIsomorphism::new(lifted_domain, u)?;
        let derived_codomain = isomorphism.codomain();

        if !ExtensionField::eq(derived_codomain.a(), lifted_twist.a())
            || !ExtensionField::eq(derived_codomain.b(), lifted_twist.b())
        {
            return Err(CurveIsomorphismError::CurvesNotIsomorphic);
        }

        Ok(isomorphism)
    }

    fn is_j_1728_family(&self) -> bool {
        F::is_zero(self.original.b())
    }

    fn certified_base_field_scaling_factor(&self) -> Option<F::Elem> {
        if let Some(square_root) = F::sqrt(&self.d) {
            return Some(square_root);
        }
        if self.is_j_1728_family() {
            let negated_factor = F::neg(&self.d);
            if let Some(square_root_of_negated_factor) = F::sqrt(&negated_factor) {
                return Some(square_root_of_negated_factor);
            }
        }
        None
    }
}

impl<F: EnumerableFiniteField + SqrtField + FiniteField> ShortWeierstrassQuadraticTwist<F> {
    /// Computes the Frobenius relation between `E` and the stored twist `E'`.
    ///
    /// If the chosen twist factor is genuinely quadratic over `F_q`, one
    /// expects `#E(F_q) + #E'(F_q) = 2q + 2`, and equivalently `t' = -t`.
    ///
    /// The current implementation derives both traces from exhaustive point
    /// counts on the two curves, records the package's base-field twist kind,
    /// and then compares the resulting invariants.
    ///
    /// Complexity: `Θ(1)`.
    pub fn frobenius_relation(&self) -> Result<QuadraticTwistFrobeniusRelation, CurveError> {
        let twist_kind = self.kind();
        let original = self.original().frobenius_trace()?;
        let twist = self.twist().frobenius_trace()?;

        let sum_orders = original.curve_order() + twist.curve_order();
        let field_order = original.field_order();
        let expected_sum = field_order * 2u8 + 2u8;
        let holds = sum_orders == expected_sum;

        Ok(QuadraticTwistFrobeniusRelation::new(
            twist_kind,
            original,
            twist,
            sum_orders,
            expected_sum,
            holds,
        ))
    }
}
