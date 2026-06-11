use std::hash::Hash;

use crate::elliptic_curves::{ShortWeierstrassCurve, ShortWeierstrassFunctionField};
use crate::fields::{EnumerableFiniteField, SqrtField};
use crate::isogenies::scalar_multiplication::ScalarMultiplicationIsogeny;
use crate::isogenies::{
    AbsoluteFrobeniusIsogeny, DualIsogenyError, FrobeniusLikeIsogeny,
    FrobeniusVerschiebungFactorizationReport, Isogeny, IsogenyConstructionError, IsogenyError,
    IsogenyMapError, ShortWeierstrassFunctionFieldMap, VerschiebungCertificate, VerschiebungError,
    VerschiebungIsogeny,
};

impl<F> ScalarMultiplicationIsogeny<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Returns the pullback map `[n]^* : F(E) -> F(E)` induced by scalar
    /// multiplication on the generic point.
    ///
    /// Let `P_gen = (x, y)` be the generic point of `E` viewed as a point of
    /// `E(F(E))`. Then the multiplication-by-`n` map is determined by the
    /// coordinates of `[n]P_gen = (X_n, Y_n)`.
    ///
    /// This method computes that generic multiple inside the existing
    /// short-Weierstrass function-field layer and returns the pullback
    /// `[n]^*(x) = X_n`, `[n]^*(y) = Y_n`.
    ///
    /// Since the constructor of [`ScalarMultiplicationIsogeny`] rejects the
    /// zero scalar, the image of the generic point is expected to stay affine
    /// in the current short-Weierstrass presentation.
    ///
    /// Complexity: `Θ(log n)` generic-point additions/doublings in `E(F(E))`
    /// from the double-and-add ladder, plus one final pullback-map
    /// validation.
    pub fn as_function_field_map(
        &self,
    ) -> Result<ShortWeierstrassFunctionFieldMap<F>, IsogenyError> {
        let field = ShortWeierstrassFunctionField::<F>::new(self.curve.clone());
        let image = field.generic_point_multiple(self.scalar())?;

        ShortWeierstrassFunctionFieldMap::new(
            self.curve.clone(),
            self.curve.clone(),
            image.x().unwrap().clone(),
            image.y().unwrap().clone(),
        )
    }

    /// Returns a certified pullback map for `[p]^*` using a verified
    /// Verschiebung certificate.
    ///
    /// This educational surface is intentionally narrower than
    /// [`Self::as_function_field_map`]:
    ///
    /// - it supports only the case `scalar = p = char(F)`
    /// - it reuses a certified factorization `[p] = V ∘ Frob_p`
    /// - it returns the corresponding pullback
    ///   `[p]^* = Frob_p^* ∘ V^*`
    ///
    /// It remains useful as an independent certification route for the
    /// characteristic-`p` map, even now that the crate can derive `[n]^*`
    /// directly from the generic point.
    ///
    /// Error policy:
    ///
    /// - returns [`IsogenyError::Dual(DualIsogenyError::DegreeMismatch)`] when `self.scalar() != p`
    /// - returns [`IsogenyError::Map(IsogenyMapError::CompositionDomainCodomainMismatch)`] when the
    ///   scalar-multiplication curve does not match the certificate's source curve `E`
    /// - returns the certificate's own duality error when the supplied
    ///   certificate is internally inconsistent
    pub fn as_function_field_map_from_verschiebung(
        &self,
        certificate: &VerschiebungCertificate<F>,
    ) -> Result<ShortWeierstrassFunctionFieldMap<F>, IsogenyError> {
        if self.scalar() != F::characteristic() {
            return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
        }

        let curve = self.domain();
        if curve != certificate.frobenius().domain()
            || curve != certificate.verschiebung().codomain_curve()
        {
            return Err(IsogenyError::Map(
                IsogenyMapError::CompositionDomainCodomainMismatch,
            ));
        }

        certificate.verify_duality_relations()?;

        let derived = certificate
            .frobenius()
            .as_function_field_map()
            .compose(certificate.verschiebung().as_function_field_map())?;

        if derived == *certificate.multiplication_by_p_on_e() {
            Ok(derived)
        } else {
            Err(IsogenyError::Verschiebung(
                VerschiebungError::LeftDualityViolation,
            ))
        }
    }

    /// Constructs the function-field-side Verschiebung `V : E^(p) -> E`
    /// directly from the pullback of `[p]`.
    ///
    /// Mathematically, in characteristic `p` one has the factorization
    ///
    /// `[p] = V \circ Frob_p`,
    ///
    /// where `Frob_p : E -> E^(p)` is absolute Frobenius and
    /// `V : E^(p) -> E` is Verschiebung. Passing to function fields reverses
    /// arrows:
    ///
    /// `[p]^* = Frob_p^* \circ V^*`.
    ///
    /// On generators of `F(E^(p))`, the Frobenius pullback acts as taking
    /// `p`-th powers, so if
    ///
    /// `V^*(x) = X_V` and `V^*(y) = Y_V`,
    ///
    /// then
    ///
    /// `[p]^*(x) = X_V^p`,
    /// `[p]^*(y) = Y_V^p`.
    ///
    /// This method first computes `[p]^*` directly from the generic point,
    /// then inverts the short-Weierstrass absolute-Frobenius pullback on its
    /// two coordinate functions, and finally interprets the recovered
    /// preimages as elements of `F(E^(p))`.
    ///
    /// Current scope note:
    ///
    /// - this reconstruction covers the absolute-Frobenius factorization
    ///   `[p] = V \circ Frob_p`
    /// - it does **not** currently cover an analogous reconstruction through
    ///   the relative Frobenius `π_q`
    ///
    /// The resulting map is then checked against the source and target curves
    /// through [`VerschiebungIsogeny::new`].
    ///
    /// Complexity: dominated by one direct construction of `[p]^*` via
    /// [`Self::as_function_field_map`], two inversions of the absolute
    /// Frobenius pullback in the function field, and one final pullback-map
    /// validation.
    pub fn verschiebung_isogeny_from_direct_p_pullback(
        &self,
    ) -> Result<VerschiebungIsogeny<F>, IsogenyError> {
        if self.scalar() != F::characteristic() {
            return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
        }

        let frobenius = AbsoluteFrobeniusIsogeny::new(self.curve.clone())?;
        let p_pullback = self.as_function_field_map()?;
        let twist_curve = frobenius.codomain().clone();

        let x_pullback = p_pullback
            .x_pullback()
            .inverse_absolute_frobenius_pullback_to_twist(&self.curve, &twist_curve)
            .ok_or(IsogenyError::Construction(
                IsogenyConstructionError::MissingInverseFrobeniusPreimageForVerschiebung {
                    coordinate: "x",
                },
            ))?;
        let y_pullback = p_pullback
            .y_pullback()
            .inverse_absolute_frobenius_pullback_to_twist(&self.curve, &twist_curve)
            .ok_or(IsogenyError::Construction(
                IsogenyConstructionError::MissingInverseFrobeniusPreimageForVerschiebung {
                    coordinate: "y",
                },
            ))?;

        let verschiebung_pullback = ShortWeierstrassFunctionFieldMap::new(
            twist_curve,
            self.curve.clone(),
            x_pullback,
            y_pullback,
        )?;

        VerschiebungIsogeny::new(frobenius, verschiebung_pullback)
    }

    /// Builds a fully verified Verschiebung certificate directly from `[p]^*`.
    ///
    /// This combines the direct generic-point pullback of `[p]` with the
    /// previous method that reconstructs `V^*` by inverting Frobenius
    /// pullbacks, and then certifies
    /// both characteristic-`p` identities against direct scalar-multiplication
    /// pullbacks on `E` and on the Frobenius twist `E^(p)`.
    ///
    /// The current implementation certifies only the absolute-Frobenius story,
    /// not a relative-Frobenius analogue.
    ///
    /// Complexity: dominated by two direct pullback constructions for `[p]`
    /// (one on `E`, one on `E^(p)`), one reconstruction of Verschiebung by
    /// inverting Frobenius pullbacks, and two pullback compositions for the
    /// certificate checks.
    pub fn verschiebung_certificate_from_direct_p_pullback(
        &self,
    ) -> Result<VerschiebungCertificate<F>, IsogenyError> {
        if self.scalar() != F::characteristic() {
            return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
        }

        let verschiebung = self.verschiebung_isogeny_from_direct_p_pullback()?;
        let multiplication_by_p_on_e = self.as_function_field_map()?;
        let multiplication_on_twist =
            ScalarMultiplicationIsogeny::new(verschiebung.domain_curve().clone(), self.scalar())?;
        let multiplication_by_p_on_frobenius_twist =
            multiplication_on_twist.as_function_field_map()?;

        VerschiebungCertificate::new(
            verschiebung,
            multiplication_by_p_on_e,
            multiplication_by_p_on_frobenius_twist,
        )
    }

    /// Packages the direct characteristic-`p` story
    ///
    /// `[p] = V \circ Frob_p`
    ///
    /// into one educational report.
    ///
    /// This report contains:
    ///
    /// - the source curve `E`
    /// - the direct pullback `[p]^*`
    /// - the absolute Frobenius `Frob_p`
    /// - the reconstructed Verschiebung `V`
    /// - the certificate of the two duality relations
    ///
    /// It is intended as the natural input to visualization helpers and
    /// examples that want to present the full factorization story from a single
    /// curve.
    ///
    /// This report is specifically about the absolute-Frobenius factorization
    /// `[p] = V \circ Frob_p`; the current crate does not yet provide the
    /// corresponding relative-Frobenius reconstruction.
    ///
    /// Complexity: the same asymptotic cost as
    /// [`Self::verschiebung_certificate_from_direct_p_pullback`] plus one
    /// direct computation of `[p]^*` on `E`.
    pub fn frobenius_verschiebung_factorization_report(
        &self,
    ) -> Result<FrobeniusVerschiebungFactorizationReport<F>, IsogenyError> {
        if self.scalar() != F::characteristic() {
            return Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch));
        }

        Ok(FrobeniusVerschiebungFactorizationReport::new(
            self.as_function_field_map()?,
            self.verschiebung_certificate_from_direct_p_pullback()?,
        ))
    }
}
