use std::hash::Hash;

use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{FiniteGroupCurveModel, GroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

use crate::isogenies::{
    DualIsogenyError, FrobeniusLikeIsogeny, Isogeny, IsogenyConstructionError, IsogenyError,
    IsogenyMapError, ShortWeierstrassFunctionFieldMap, VerschiebungCertificate, VerschiebungError,
};

/// Scalar-multiplication isogeny `[n] : E -> E` on a small explicit curve group.
///
/// For a non-zero integer `n`, elliptic-curve multiplication by `n`
///
/// `[n](P) = P + P + ... + P`
///
/// is a standard example of an isogeny from a curve to itself.
///
/// In this educational implementation:
///
/// - the domain and codomain are the same curve value
/// - the degree is reported as `n^2`
/// - `kernel_points()` exposes the points of `E(F_q)` killed by `[n]`
pub struct ScalarMultiplicationIsogeny<C: GroupCurveModel> {
    curve: C,
    scalar: u64,
    kernel_points: Vec<C::Point>,
}

impl<C: FiniteGroupCurveModel> ScalarMultiplicationIsogeny<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    /// Builds the scalar-multiplication isogeny `[n]`.
    ///
    /// The current constructor is intentionally restricted to small finite
    /// curve groups so it can materialize the kernel
    ///
    /// `E[n] = { P in E(F_q) : [n]P = O }`.
    ///
    /// Scalar `0` is rejected because this crate reserves the isogeny surface
    /// for the usual non-constant multiplication-by-`n` maps.
    pub fn new(curve: C, scalar: u64) -> Result<Self, IsogenyError> {
        if scalar == 0 {
            return Err(IsogenyError::Construction(
                IsogenyConstructionError::ZeroScalarIsNotIsogeny,
            ));
        }

        let identity = curve.identity();
        let kernel_points = curve
            .points()
            .into_iter()
            .map(|point| -> Result<_, IsogenyError> {
                let image = curve.mul_scalar(&point, scalar)?;
                Ok((point, image == identity))
            })
            .collect::<Result<Vec<_>, IsogenyError>>()?
            .into_iter()
            .filter_map(|(point, kills_point)| kills_point.then_some(point))
            .collect();

        Ok(Self {
            curve,
            scalar,
            kernel_points,
        })
    }

    /// Returns the underlying scalar `n` in `[n]`.
    pub fn scalar(&self) -> u64 {
        self.scalar
    }
}

impl<C: GroupCurveModel> Isogeny<C, C> for ScalarMultiplicationIsogeny<C>
where
    C::Point: Clone,
{
    fn domain(&self) -> &C {
        &self.curve
    }

    fn codomain(&self) -> &C {
        &self.curve
    }

    fn degree(&self) -> usize {
        usize::try_from(u128::from(self.scalar) * u128::from(self.scalar))
            .expect("educational scalar-multiplication degrees should fit in usize")
    }

    fn evaluate(&self, p: &C::Point) -> Result<C::Point, IsogenyError> {
        self.curve.mul_scalar(p, self.scalar).map_err(Into::into)
    }

    fn kernel_points(&self) -> &[C::Point] {
        &self.kernel_points
    }
}

impl<F> ScalarMultiplicationIsogeny<ShortWeierstrassCurve<F>>
where
    F: EnumerableFiniteField + SqrtField + Clone,
    F::Elem: Clone + Eq + Hash + PartialEq,
{
    /// Returns a certified pullback map for `[p]^*` using a verified
    /// Verschiebung certificate.
    ///
    /// This current educational surface is intentionally narrow:
    ///
    /// - it supports only the case `scalar = p = char(F)`
    /// - it reuses a certified factorization `[p] = V ∘ Frob_p`
    /// - it returns the corresponding pullback
    ///   `[p]^* = Frob_p^* ∘ V^*`
    ///
    /// It does **not** yet derive the rational pullback of `[n]` for arbitrary
    /// scalars from point evaluation alone.
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
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::elliptic_curves::{
        AffineCurveModel, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel, GroupCurveModel,
        ShortWeierstrassCurve,
    };
    use crate::fields::{Field, Fp};
    use crate::isogenies::{
        AbsoluteFrobeniusIsogeny, DualIsogenyError, FrobeniusLikeIsogeny, Isogeny,
        IsogenyConstructionError, IsogenyError, ScalarMultiplicationIsogeny, VerifiableIsogeny,
        VerschiebungCertificate, VerschiebungIsogeny,
    };

    type F41 = Fp<41>;
    type Curve = ShortWeierstrassCurve<F41>;

    fn curve() -> Curve {
        Curve::new(F41::from_i64(2), F41::from_i64(3)).expect("valid curve")
    }

    #[test]
    fn constructor_rejects_zero_scalar() {
        assert!(matches!(
            ScalarMultiplicationIsogeny::new(curve(), 0),
            Err(IsogenyError::Construction(
                IsogenyConstructionError::ZeroScalarIsNotIsogeny
            ))
        ));
    }

    #[test]
    fn degree_of_multiplication_by_two_is_four() {
        let isogeny =
            ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

        assert_eq!(isogeny.degree(), 4);
        assert_eq!(isogeny.scalar(), 2);
    }

    #[test]
    fn evaluation_matches_group_scalar_multiplication() {
        let curve = curve();
        let point = curve
            .point(F41::from_i64(3), F41::from_i64(6))
            .expect("sample point should lie on the curve");
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 3)
            .expect("scalar isogeny should build");

        assert_eq!(
            isogeny.evaluate(&point),
            curve.mul_scalar(&point, 3).map_err(Into::into)
        );
    }

    #[test]
    fn scalar_one_is_identity_map() {
        let curve = curve();
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 1)
            .expect("scalar-one isogeny should build");

        for point in curve.points() {
            assert_eq!(
                isogeny
                    .evaluate(&point)
                    .expect("scalar-one isogeny should evaluate"),
                point
            );
        }
    }

    #[test]
    fn kernel_points_match_the_rational_two_torsion_plus_identity() {
        let curve = curve();
        let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), 2)
            .expect("scalar isogeny should build");

        let mut expected = vec![curve.identity()];
        expected.extend(curve.points_of_order(2));

        assert_eq!(isogeny.kernel_points(), expected.as_slice());
    }

    #[test]
    fn exhaustive_verifier_passes_for_multiplication_by_two() {
        let isogeny =
            ScalarMultiplicationIsogeny::new(curve(), 2).expect("scalar isogeny should build");

        assert_eq!(isogeny.verify_maps_domain_to_codomain(), Ok(()));
        assert_eq!(isogeny.verify_maps_kernel_to_identity(), Ok(()));
        assert_eq!(isogeny.verify_homomorphism(), Ok(()));
        assert_eq!(isogeny.verify_kernel_exactness(), Ok(()));
    }

    #[test]
    fn function_field_map_from_verschiebung_recovers_the_certified_p_pullback() {
        let curve = curve();
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
        let certificate =
            VerschiebungCertificate::new(verschiebung, expected_left.clone(), expected_right)
                .expect("certificate should build");

        let scalar = ScalarMultiplicationIsogeny::new(curve, 41)
            .expect("scalar multiplication should build");

        assert_eq!(
            scalar
                .as_function_field_map_from_verschiebung(&certificate)
                .expect("certified map should build"),
            expected_left
        );
    }

    #[test]
    fn function_field_map_from_verschiebung_rejects_non_characteristic_scalar() {
        let curve = curve();
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

        let scalar =
            ScalarMultiplicationIsogeny::new(curve, 2).expect("scalar multiplication should build");

        assert_eq!(
            scalar.as_function_field_map_from_verschiebung(&certificate),
            Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch))
        );
    }

    fn curve_and_point() -> impl Strategy<Value = (Curve, <Curve as CurveModel>::Point)> {
        let curve = curve();
        let points = curve.points();
        let len = points.len();

        (0usize..len).prop_map(move |index| (curve.clone(), points[index].clone()))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(24))]

        #[test]
        fn property_scalar_isogeny_evaluation_matches_curve_scalar_multiplication(
            (curve, point) in curve_and_point(),
            scalar in 1u64..6,
        ) {
            let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), scalar)
                .expect("scalar isogeny should build");

            prop_assert_eq!(
                isogeny.evaluate(&point),
                curve.mul_scalar(&point, scalar).map_err(Into::into)
            );
        }

        #[test]
        fn property_scalar_isogeny_kernel_matches_points_killed_by_the_scalar(
            scalar in 1u64..6,
        ) {
            let curve = curve();
            let isogeny = ScalarMultiplicationIsogeny::new(curve.clone(), scalar)
                .expect("scalar isogeny should build");
            let expected: Vec<_> = curve
                .points()
                .into_iter()
                .filter(|point| curve.mul_scalar(point, scalar).ok() == Some(curve.identity()))
                .collect();

            prop_assert_eq!(isogeny.kernel_points(), expected.as_slice());
        }
    }
}
