use super::shared::curve;
use crate::elliptic_curves::short_weierstrass::isogenies::{
    frobenius::{AbsoluteFrobeniusIsogeny, VerschiebungCertificate, VerschiebungIsogeny},
    function_field_maps::ShortWeierstrassFunctionFieldMap,
};
use crate::elliptic_curves::{
    ShortWeierstrassCurve, short_weierstrass::isogenies::frobenius::FrobeniusLikeIsogeny,
};
use crate::fields::extension_field::define_fp_quadratic_extension;
use crate::fields::traits::*;
use crate::isogenies::{
    error::{DualIsogenyError, IsogenyError},
    scalar_multiplication::ScalarMultiplicationIsogeny,
    traits::Isogeny,
};

define_fp_quadratic_extension!(
    spec: F5Sqrt2ScalarMultiplicationSpec,
    field: F5Sqrt2ScalarMultiplication,
    base: crate::fields::Fp5,
    non_residue: 2,
    name: "F5(sqrt(2)) for scalar-multiplication Frobenius tests",
);

fn nontrivial_extension_curve() -> ShortWeierstrassCurve<F5Sqrt2ScalarMultiplication> {
    let alpha = F5Sqrt2ScalarMultiplication::element(vec![
        crate::fields::Fp5::zero(),
        crate::fields::Fp5::one(),
    ]);
    ShortWeierstrassCurve::<F5Sqrt2ScalarMultiplication>::new(
        alpha,
        F5Sqrt2ScalarMultiplication::one(),
    )
    .expect("valid curve over F5^2")
}

fn multiplication_by_p_on_curve()
-> ScalarMultiplicationIsogeny<ShortWeierstrassCurve<crate::fields::Fp41>> {
    ScalarMultiplicationIsogeny::new(curve(), 41).expect("scalar multiplication should build")
}

fn certified_verschiebung_fixture(
    curve: &ShortWeierstrassCurve<crate::fields::Fp41>,
) -> (
    VerschiebungCertificate<crate::fields::Fp41>,
    ShortWeierstrassFunctionFieldMap<crate::fields::Fp41>,
) {
    let frobenius =
        AbsoluteFrobeniusIsogeny::new(curve.clone()).expect("absolute Frobenius should build");
    let candidate_v = ShortWeierstrassFunctionFieldMap::new(
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
    let verschiebung =
        VerschiebungIsogeny::new(frobenius.clone(), candidate_v).expect("candidate should build");
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

    (certificate, expected_left)
}

#[test]
fn function_field_map_from_verschiebung_recovers_the_certified_p_pullback() {
    let curve = curve();
    let (certificate, expected_left) = certified_verschiebung_fixture(&curve);
    let scalar = multiplication_by_p_on_curve();

    assert_eq!(
        scalar
            .as_function_field_map_from_verschiebung(&certificate)
            .expect("certified map should build"),
        expected_left
    );
}

#[test]
fn direct_p_pullback_can_build_a_verschiebung_isogeny() {
    let scalar = multiplication_by_p_on_curve();
    let verschiebung = scalar
        .verschiebung_isogeny_from_direct_p_pullback()
        .expect("Verschiebung should be extracted from [p]^*");

    assert_eq!(verschiebung.codomain_curve(), scalar.domain());
    assert_eq!(
        verschiebung.domain_curve(),
        verschiebung.frobenius().codomain()
    );
    assert_eq!(verschiebung.degree(), 41);
}

#[test]
fn direct_p_pullback_can_build_a_verified_verschiebung_over_nontrivial_extension_curve() {
    let curve = nontrivial_extension_curve();
    let scalar = ScalarMultiplicationIsogeny::new(curve.clone(), 5)
        .expect("scalar multiplication should build");
    let verschiebung = scalar
        .verschiebung_isogeny_from_direct_p_pullback()
        .expect("Verschiebung should be extracted from [p]^*");
    let certificate = scalar
        .verschiebung_certificate_from_direct_p_pullback()
        .expect("certificate should build");

    assert_eq!(verschiebung.codomain_curve(), &curve);
    assert_eq!(
        verschiebung.domain_curve(),
        verschiebung.frobenius().codomain()
    );
    assert_ne!(verschiebung.domain_curve(), verschiebung.codomain_curve());
    assert_eq!(certificate.verify_duality_relations(), Ok(()));
}

#[test]
fn function_field_map_from_verschiebung_rejects_non_characteristic_scalar() {
    let curve = curve();
    let (certificate, _) = certified_verschiebung_fixture(&curve);

    let scalar =
        ScalarMultiplicationIsogeny::new(curve, 2).expect("scalar multiplication should build");

    assert_eq!(
        scalar.as_function_field_map_from_verschiebung(&certificate),
        Err(IsogenyError::Dual(DualIsogenyError::DegreeMismatch))
    );
}
