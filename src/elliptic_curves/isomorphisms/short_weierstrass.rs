use core::fmt;

use crate::{
    elliptic_curves::{
        AffinePoint, CurveIsomorphism, CurveIsomorphismError, CurveModel, ShortWeierstrassCurve,
    },
    fields::{ExtensionField, ExtensionFieldSpec, Field, FieldError, PolynomialModulus, SqrtField},
};

/// Explicit short-Weierstrass scaling isomorphism determined by a parameter
/// `u`.
///
/// The intended mathematical convention is the map
/// `\phi_u : E -> E'`
/// defined on affine points by
/// `(x, y) -> (u^2 x, u^3 y)`.
///
/// If the domain curve is written in short-Weierstrass form as
/// `E: y^2 = x^3 + ax + b`,
/// then the image curve is
/// `E': y^2 = x^3 + a'x + b'`
/// with transformed coefficients
/// `a' = u^4 a`, `b' = u^6 b`.
///
/// This type treats `domain` and the parameter `u` as the primary data. The
/// codomain is derived automatically from those values instead of being stored
/// as a second source of truth.
pub struct ShortWeierstrassIsomorphism<F: Field> {
    domain: ShortWeierstrassCurve<F>,
    codomain: ShortWeierstrassCurve<F>,
    u: F::Elem,
}

impl<F: Field> Clone for ShortWeierstrassIsomorphism<F> {
    fn clone(&self) -> Self {
        Self {
            domain: self.domain.clone(),
            codomain: self.codomain.clone(),
            u: self.u.clone(),
        }
    }
}

impl<F: Field> fmt::Debug for ShortWeierstrassIsomorphism<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ShortWeierstrassIsomorphism")
            .field("domain", &self.domain)
            .field("codomain", &self.codomain)
            .field("u", &self.u)
            .finish()
    }
}

impl<F: Field> ShortWeierstrassIsomorphism<F> {
    /// Builds the short-Weierstrass scaling isomorphism determined by `u`.
    ///
    /// Since `\phi_u` must be invertible, `u` must be invertible in the base
    /// field.
    pub fn new(
        domain: ShortWeierstrassCurve<F>,
        u: F::Elem,
    ) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&u).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }

        let codomain = Self::derive_codomain_from(&domain, &u)?;
        let isomorphism = Self {
            domain,
            codomain,
            u,
        };

        Ok(isomorphism)
    }

    /// Returns the scaling parameter `u` from `\phi_u`.
    pub fn scaling_factor(&self) -> &F::Elem {
        &self.u
    }

    /// Returns the inverse isomorphism `\phi_u^{-1}`.
    ///
    /// Under the convention
    /// `\phi_u(x, y) = (u^2 x, u^3 y)`,
    /// the inverse map is the short-Weierstrass scaling determined by `u^{-1}`
    /// from the derived codomain `E'` back to the original domain `E`.
    pub fn inverse(&self) -> Result<Self, CurveIsomorphismError> {
        let inverse_u = F::inv(&self.u).ok_or(CurveIsomorphismError::NonInvertibleScale)?;
        Self::new(self.codomain.clone(), inverse_u)
    }

    fn derive_codomain_from(
        domain: &ShortWeierstrassCurve<F>,
        u: &F::Elem,
    ) -> Result<ShortWeierstrassCurve<F>, CurveIsomorphismError> {
        let u2 = F::square(u);
        let u4 = F::square(&u2);
        let u6 = F::mul(&u4, &u2);

        ShortWeierstrassCurve::new(F::mul(&u4, domain.a()), F::mul(&u6, domain.b()))
            .map_err(Into::into)
    }

    fn u2(&self) -> F::Elem {
        F::square(&self.u)
    }

    fn u3(&self) -> F::Elem {
        F::mul(&self.u2(), &self.u)
    }
}

impl<F: Field> CurveIsomorphism for ShortWeierstrassIsomorphism<F> {
    type Domain = ShortWeierstrassCurve<F>;
    type Codomain = ShortWeierstrassCurve<F>;

    fn domain(&self) -> &Self::Domain {
        &self.domain
    }

    fn codomain(&self) -> &Self::Codomain {
        &self.codomain
    }

    fn evaluate(
        &self,
        point: &<Self::Domain as CurveModel>::Point,
    ) -> Result<<Self::Codomain as CurveModel>::Point, CurveIsomorphismError> {
        if !self.domain.contains(point) {
            return Err(CurveIsomorphismError::PointNotOnDomain);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::infinity()),
            AffinePoint::Finite { x, y } => {
                let image = AffinePoint::new(F::mul(&self.u2(), x), F::mul(&self.u3(), y));

                if !self.codomain.contains(&image) {
                    return Err(CurveIsomorphismError::ImagePointNotOnCodomain);
                }

                Ok(image)
            }
        }
    }
}

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
        PolynomialModulus::<F>::new(vec![F::neg(&self.d), F::zero(), F::one()]).map_err(Into::into)
    }

    fn lift_curve_to_extension<S>(
        curve: &ShortWeierstrassCurve<F>,
    ) -> Result<ShortWeierstrassCurve<ExtensionField<S>>, CurveIsomorphismError>
    where
        S: ExtensionFieldSpec<Base = F>,
    {
        ShortWeierstrassCurve::<ExtensionField<S>>::new(
            ExtensionField::<S>::from_base(curve.a().clone()),
            ExtensionField::<S>::from_base(curve.b().clone()),
        )
        .map_err(Into::into)
    }
}

impl<F: SqrtField> ShortWeierstrassQuadraticTwist<F> {
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

        // No extra `j = 0` branch is needed for this quadratic-twist
        // normalization. If `a = 0`, then a base-field scaling from
        // `E : y^2 = x^3 + b` to `E^(d) : y^2 = x^3 + d^3 b` would satisfy
        // `u^6 = d^3`. Since `d != 0` for a valid twist factor, this already
        // forces `d = (d^2 / u^3)^2`, so `d` must be a square in the base
        // field. The ordinary square-root witness is therefore exact even
        // though the geometric automorphism group is larger at `j = 0`.

        None
    }

    /// Returns whether the stored twist is trivial or genuinely quadratic over
    /// the current base field.
    ///
    /// Generic curves still use the usual square-root witness for the twist
    /// factor `d`.
    ///
    /// For the special `j = 1728` family `y^2 = x^3 + ax`, the implementation
    /// also certifies the extra base-field-trivial path coming from a square
    /// root of `-d`, since `u^2 = -d` implies `u^4 = d^2` and therefore
    /// identifies `E` with `E^(d)` when `b = 0`.
    ///
    /// The other exceptional family `j = 0` (`a = 0`) also has extra
    /// geometric automorphisms, but for the specific quadratic-twist
    /// normalization `b' = d^3 b` no extra base-field-trivial path appears:
    /// any base-field scaling would satisfy `u^6 = d^3`, which already forces
    /// `d` itself to be a square. So the ordinary square-root witness remains
    /// the whole story there.
    pub fn kind(&self) -> TwistKind {
        if self.certified_base_field_scaling_factor().is_some() {
            TwistKind::Trivial
        } else {
            TwistKind::Quadratic
        }
    }

    /// Returns one certified base-field isomorphism from `E` to `E^(d)` when
    /// the current field backend can witness it directly.
    ///
    /// In the generic case, this comes from a square root of `d`.
    ///
    /// For the exceptional `j = 1728` family, this helper also recognizes the
    /// extra base-field-trivial path coming from a square root of `-d`.
    ///
    /// For the exceptional `j = 0` family, no additional branch is needed for
    /// this quadratic-twist normalization: the existence of a base-field
    /// scaling already forces `d` to be a square.
    pub fn base_field_isomorphism(&self) -> Option<ShortWeierstrassIsomorphism<F>> {
        let scaling_factor = self.certified_base_field_scaling_factor()?;
        ShortWeierstrassIsomorphism::new(self.original_curve_copy(), scaling_factor).ok()
    }

    /// Returns the canonical scaling isomorphism over a genuine quadratic
    /// extension presented as `F[x] / (x^2 - d)`.
    ///
    /// This helper is intentionally restricted to the genuinely quadratic
    /// case. If the stored twist is already base-field trivial, whether by the
    /// ordinary square-root witness or by the exceptional `j = 1728`
    /// base-field path, then [`Self::base_field_isomorphism`] should be used
    /// instead.
    ///
    /// The caller must supply an [`ExtensionFieldSpec`] whose defining modulus
    /// is exactly `x^2 - d`. The resulting isomorphism uses the class of `x`
    /// in that quotient as the scaling witness `u`, so `u^2 = d` by
    /// construction.
    pub fn isomorphism_over_quadratic_extension<S>(
        &self,
    ) -> Result<ShortWeierstrassIsomorphism<ExtensionField<S>>, CurveIsomorphismError>
    where
        S: ExtensionFieldSpec<Base = F>,
    {
        if self.kind() == TwistKind::Trivial {
            return Err(FieldError::NonIrreduciblePolynomial.into());
        }

        let expected_modulus = self.expected_quadratic_extension_modulus()?;
        if S::defining_modulus() != expected_modulus {
            return Err(FieldError::IncompatibleFieldParameters.into());
        }

        ExtensionField::<S>::check_structure()?;

        let lifted_domain = Self::lift_curve_to_extension::<S>(&self.original)?;
        let lifted_twist = Self::lift_curve_to_extension::<S>(&self.twist)?;
        let u = ExtensionField::<S>::element(vec![F::zero(), F::one()]);
        let isomorphism = ShortWeierstrassIsomorphism::new(lifted_domain, u)?;
        let derived_codomain = isomorphism.codomain();

        if !ExtensionField::<S>::eq(derived_codomain.a(), lifted_twist.a())
            || !ExtensionField::<S>::eq(derived_codomain.b(), lifted_twist.b())
        {
            return Err(CurveIsomorphismError::CurvesNotIsomorphic);
        }

        Ok(isomorphism)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::elliptic_curves::isomorphisms::{
        CurveIsomorphism, ShortWeierstrassIsomorphism, ShortWeierstrassQuadraticTwist, TwistKind,
    };
    use crate::proptest_support::config::CurveStrategyConfig;
    use crate::proptest_support::isogenies::arb_short_weierstrass_isomorphism_case;
    use crate::{
        elliptic_curves::{
            AffineCurveModel, AffinePoint, CurveIsomorphismError, CurveModel, ShortWeierstrassCurve,
        },
        fields::{CbrtField, EnumerableFiniteField, Field, FieldError, Fp, SqrtField},
    };

    type F7 = Fp<7>;
    type F13 = Fp<13>;
    type F19 = Fp<19>;

    crate::fields::define_fp_quadratic_extension!(
        spec: F7Sqrt3Spec,
        field: F7Sqrt3,
        base: F7,
        non_residue: 3,
        name: "F7(sqrt(3))",
    );

    #[allow(dead_code)]
    type _F7Sqrt3Marker = F7Sqrt3;

    crate::fields::define_fp_quadratic_extension!(
        spec: F19Sqrt2Spec,
        field: F19Sqrt2,
        base: F19,
        non_residue: 2,
        name: "F19(sqrt(2))",
    );

    crate::fields::define_fp_quadratic_extension!(
        spec: F19Sqrt3Spec,
        field: F19Sqrt3,
        base: F19,
        non_residue: 3,
        name: "F19(sqrt(3))",
    );

    fn f7_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3)).expect("valid curve")
    }

    fn f19_curve() -> ShortWeierstrassCurve<F19> {
        ShortWeierstrassCurve::<F19>::new(F19::from_i64(2), F19::from_i64(3)).expect("valid curve")
    }

    fn f7_j1728_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(1), F7::zero()).expect("valid j=1728 curve")
    }

    fn f13_j0_curve() -> ShortWeierstrassCurve<F13> {
        ShortWeierstrassCurve::<F13>::new(F13::zero(), F13::from_i64(1)).expect("valid j=0 curve")
    }

    fn f7_point(x: i64, y: i64) -> AffinePoint<F7> {
        f7_curve()
            .point(F7::from_i64(x), F7::from_i64(y))
            .expect("valid point on the sample curve")
    }

    fn f19_point(x: i64, y: i64) -> AffinePoint<F19> {
        f19_curve()
            .point(F19::from_i64(x), F19::from_i64(y))
            .expect("valid point on the sample curve")
    }

    #[test]
    fn constructor_rejects_noninvertible_scale() {
        assert!(matches!(
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::zero()),
            Err(CurveIsomorphismError::NonInvertibleScale)
        ));
    }

    #[test]
    fn domain_getter_returns_original_curve() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

        assert!(F7::eq(isomorphism.domain().a(), &F7::from_i64(2)));
        assert!(F7::eq(isomorphism.domain().b(), &F7::from_i64(3)));
    }

    #[test]
    fn codomain_is_derived_from_u4_and_u6() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let codomain = isomorphism.codomain();

        assert!(F7::eq(codomain.a(), &F7::from_i64(1)));
        assert!(F7::eq(codomain.b(), &F7::from_i64(3)));
    }

    #[test]
    fn scaling_factor_getter_returns_u() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

        assert!(F7::eq(isomorphism.scaling_factor(), &F7::from_i64(3)));
    }

    #[test]
    fn inverse_uses_inverse_scaling_factor() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let inverse = isomorphism.inverse().expect("inverse should exist");

        assert!(F7::eq(inverse.scaling_factor(), &F7::from_i64(5)));
    }

    #[test]
    fn inverse_domain_is_the_original_codomain() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let codomain = isomorphism.codomain();
        let inverse = isomorphism.inverse().expect("inverse should exist");

        assert!(F7::eq(inverse.domain().a(), codomain.a()));
        assert!(F7::eq(inverse.domain().b(), codomain.b()));
    }

    #[test]
    fn inverse_codomain_is_the_original_domain() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let inverse = isomorphism.inverse().expect("inverse should exist");
        let inverse_codomain = inverse.codomain();

        assert!(F7::eq(inverse_codomain.a(), isomorphism.domain().a()));
        assert!(F7::eq(inverse_codomain.b(), isomorphism.domain().b()));
    }

    #[test]
    fn inverse_of_inverse_returns_original_scaling_and_curves() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let double_inverse = isomorphism
            .inverse()
            .expect("inverse should exist")
            .inverse()
            .expect("double inverse should exist");

        assert!(F7::eq(
            double_inverse.scaling_factor(),
            isomorphism.scaling_factor()
        ));
        assert!(F7::eq(
            double_inverse.domain().a(),
            isomorphism.domain().a()
        ));
        assert!(F7::eq(
            double_inverse.domain().b(),
            isomorphism.domain().b()
        ));

        let double_inverse_codomain = double_inverse.codomain();
        let original_codomain = isomorphism.codomain();
        assert!(F7::eq(double_inverse_codomain.a(), original_codomain.a()));
        assert!(F7::eq(double_inverse_codomain.b(), original_codomain.b()));
    }

    #[test]
    fn evaluate_sends_infinity_to_infinity() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");

        assert_eq!(
            isomorphism
                .evaluate(&AffinePoint::infinity())
                .expect("the identity should map to itself"),
            AffinePoint::infinity()
        );
    }

    #[test]
    fn evaluate_rejects_point_outside_the_domain() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let invalid = AffinePoint::<F7>::new(F7::from_i64(2), F7::from_i64(2));

        assert!(matches!(
            isomorphism.evaluate(&invalid),
            Err(CurveIsomorphismError::PointNotOnDomain)
        ));
    }

    #[test]
    fn evaluate_transports_finite_points_to_the_codomain() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let point = f7_point(2, 1);
        let image = isomorphism
            .evaluate(&point)
            .expect("a domain point should map into the codomain");

        assert!(isomorphism.codomain().contains(&image));
        assert_eq!(image, AffinePoint::new(F7::from_i64(4), F7::from_i64(6)));
    }

    #[test]
    fn inverse_recovers_the_original_point_after_evaluation() {
        let isomorphism =
            ShortWeierstrassIsomorphism::new(f7_curve(), F7::from_i64(3)).expect("valid scaling");
        let point = f7_point(2, 1);
        let image = isomorphism
            .evaluate(&point)
            .expect("a domain point should map into the codomain");

        assert_eq!(
            isomorphism
                .inverse()
                .expect("inverse should exist")
                .evaluate(&image)
                .expect("the inverse should recover the original point"),
            point
        );
    }

    #[test]
    fn quadratic_twist_package_stores_original_twist_and_factor() {
        let original = f19_curve();
        let package = ShortWeierstrassQuadraticTwist::new(original, F19::from_i64(2))
            .expect("non-zero twist factor should produce a valid package");

        assert!(F19::eq(package.original().a(), &F19::from_i64(2)));
        assert!(F19::eq(package.original().b(), &F19::from_i64(3)));
        assert!(F19::eq(package.factor(), &F19::from_i64(2)));
        assert!(package.original().has_same_j_invariant(package.twist()));
    }

    #[test]
    fn quadratic_twist_kind_is_trivial_when_the_factor_is_a_square() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
            .expect("square twist factor should produce a valid package");

        assert_eq!(package.kind(), TwistKind::Trivial);
    }

    #[test]
    fn quadratic_twist_kind_is_quadratic_when_the_factor_is_not_a_square() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
            .expect("non-square twist factor should produce a valid package");

        assert_eq!(package.kind(), TwistKind::Quadratic);
    }

    #[test]
    fn base_field_isomorphism_exists_for_generic_square_and_non_square_factors() {
        let trivial_package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
            .expect("square twist factor should produce a valid package");
        let quadratic_package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
            .expect("non-square twist factor should produce a valid package");

        assert!(trivial_package.base_field_isomorphism().is_some());
        assert!(quadratic_package.base_field_isomorphism().is_none());
    }

    #[test]
    fn j1728_non_square_factor_can_still_be_trivial() {
        let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
            .expect("non-square twist factor should still produce a valid package");

        assert_eq!(package.kind(), TwistKind::Trivial);
    }

    #[test]
    fn j1728_non_square_factor_can_still_produce_a_base_field_isomorphism() {
        let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
            .expect("non-square twist factor should still produce a valid package");
        let isomorphism = package
            .base_field_isomorphism()
            .expect("j = 1728 should admit the extra base-field witness here");

        assert!(F7::eq(isomorphism.scaling_factor(), &F7::from_i64(2)));
        assert!(F7::eq(isomorphism.codomain().a(), package.twist().a()));
        assert!(F7::eq(isomorphism.codomain().b(), package.twist().b()));
    }

    #[test]
    fn quadratic_extension_isomorphism_rejects_j1728_twists_that_are_already_trivial() {
        let package = ShortWeierstrassQuadraticTwist::new(f7_j1728_curve(), F7::from_i64(3))
            .expect("non-square twist factor should still produce a valid package");

        assert!(matches!(
            package.isomorphism_over_quadratic_extension::<F7Sqrt3Spec>(),
            Err(CurveIsomorphismError::Field(
                FieldError::NonIrreduciblePolynomial
            ))
        ));
    }

    #[test]
    fn j0_square_factor_is_still_trivial() {
        let package = ShortWeierstrassQuadraticTwist::new(f13_j0_curve(), F13::from_i64(4))
            .expect("square twist factor should produce a valid package");

        assert_eq!(package.kind(), TwistKind::Trivial);
        assert!(package.base_field_isomorphism().is_some());
    }

    #[test]
    fn j0_non_square_factor_stays_quadratic_in_the_sample_prime_field() {
        let package = ShortWeierstrassQuadraticTwist::new(f13_j0_curve(), F13::from_i64(2))
            .expect("non-square twist factor should produce a valid package");

        assert_eq!(package.kind(), TwistKind::Quadratic);
        assert!(package.base_field_isomorphism().is_none());
    }

    #[test]
    fn j0_sample_field_has_nontrivial_cube_roots_of_unity_but_they_are_still_squares() {
        let nontrivial_cube_roots_of_unity = F13::elements()
            .into_iter()
            .filter(|element| {
                F13::eq(&F13::cube(element), &F13::one()) && !F13::eq(element, &F13::one())
            })
            .collect::<Vec<_>>();

        assert!(!nontrivial_cube_roots_of_unity.is_empty());
        assert!(
            nontrivial_cube_roots_of_unity
                .iter()
                .all(F13::has_square_root)
        );
        assert!(
            nontrivial_cube_roots_of_unity
                .iter()
                .all(|element| !F13::has_cube_root(element))
        );
    }

    #[test]
    fn base_field_isomorphism_transports_points_into_the_stored_twist() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
            .expect("square twist factor should produce a valid package");
        let point = f19_point(3, 6);
        let isomorphism = package
            .base_field_isomorphism()
            .expect("a square twist factor should produce a base-field isomorphism");
        let image = isomorphism
            .evaluate(&point)
            .expect("the base-field isomorphism should transport points");

        assert!(package.twist().contains(&image));
        assert_eq!(
            isomorphism
                .inverse()
                .expect("inverse should exist")
                .evaluate(&image)
                .expect("inverse should recover the original point"),
            point
        );
    }

    #[test]
    fn quadratic_extension_isomorphism_uses_the_expected_extension_and_codomain() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
            .expect("non-square twist factor should produce a valid package");
        let isomorphism = package
            .isomorphism_over_quadratic_extension::<F19Sqrt2Spec>()
            .expect("the matching quadratic extension should produce an isomorphism");
        let lifted_twist = ShortWeierstrassCurve::<F19Sqrt2>::new(
            F19Sqrt2::from_base(*package.twist().a()),
            F19Sqrt2::from_base(*package.twist().b()),
        )
        .expect("the stored twist should lift to the extension field");
        let derived_codomain = isomorphism.codomain();

        assert!(F19Sqrt2::eq(
            &F19Sqrt2::square(isomorphism.scaling_factor()),
            &F19Sqrt2::from_base(F19::from_i64(2))
        ));
        assert!(F19Sqrt2::eq(derived_codomain.a(), lifted_twist.a()));
        assert!(F19Sqrt2::eq(derived_codomain.b(), lifted_twist.b()));
    }

    #[test]
    fn quadratic_extension_isomorphism_transports_points_and_inverse_recovers_them() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
            .expect("non-square twist factor should produce a valid package");
        let isomorphism = package
            .isomorphism_over_quadratic_extension::<F19Sqrt2Spec>()
            .expect("the matching quadratic extension should produce an isomorphism");
        let lifted_domain = ShortWeierstrassCurve::<F19Sqrt2>::new(
            F19Sqrt2::from_base(*package.original().a()),
            F19Sqrt2::from_base(*package.original().b()),
        )
        .expect("the original curve should lift to the extension field");
        let point = lifted_domain
            .point(
                F19Sqrt2::from_base(F19::from_i64(3)),
                F19Sqrt2::from_base(F19::from_i64(6)),
            )
            .expect("the sample point should stay on the lifted curve");
        let image = isomorphism
            .evaluate(&point)
            .expect("the quadratic-extension isomorphism should transport points");

        assert!(isomorphism.codomain().contains(&image));
        assert_eq!(
            isomorphism
                .inverse()
                .expect("inverse should exist")
                .evaluate(&image)
                .expect("inverse should recover the original point"),
            point
        );
    }

    #[test]
    fn quadratic_extension_isomorphism_rejects_square_twist_factors() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(4))
            .expect("square twist factor should produce a valid package");

        assert!(matches!(
            package.isomorphism_over_quadratic_extension::<F19Sqrt2Spec>(),
            Err(CurveIsomorphismError::Field(
                FieldError::NonIrreduciblePolynomial
            ))
        ));
    }

    #[test]
    fn quadratic_extension_isomorphism_rejects_incompatible_extension_specs() {
        let package = ShortWeierstrassQuadraticTwist::new(f19_curve(), F19::from_i64(2))
            .expect("non-square twist factor should produce a valid package");
        let _field = F19Sqrt3::new().expect("the alternative quadratic extension should validate");

        assert!(matches!(
            package.isomorphism_over_quadratic_extension::<F19Sqrt3Spec>(),
            Err(CurveIsomorphismError::Field(
                FieldError::IncompatibleFieldParameters
            ))
        ));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(20))]

        #[test]
        fn generated_short_weierstrass_isomorphisms_preserve_domain_points(
            case in arb_short_weierstrass_isomorphism_case::<19>(CurveStrategyConfig::default()),
        ) {
            let image = case
                .isomorphism
                .evaluate(&case.sample_point)
                .expect("generated domain point should evaluate");

            prop_assert!(case.curve.contains(&case.sample_point));
            prop_assert!(case.isomorphism.codomain().contains(&image));
            prop_assert_eq!(
                case.isomorphism
                    .inverse()
                    .expect("generated isomorphism should be invertible")
                    .evaluate(&image)
                    .expect("inverse should recover the sampled point"),
                case.sample_point
            );
        }
    }
}
