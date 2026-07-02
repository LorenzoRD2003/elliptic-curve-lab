use super::shared::curve;
use crate::elliptic_curves::traits::{CurveModel, FiniteGroupCurveModel};
use crate::isogenies::{
    kernel::KernelDescription, scalar_multiplication::ScalarMultiplicationIsogeny, traits::Isogeny,
};

#[test]
fn characteristic_divisible_scalar_reports_mixed_kernel_description() {
    let curve = curve();
    let isogeny = ScalarMultiplicationIsogeny::new(curve, 41).expect("scalar isogeny should build");

    let description = isogeny.kernel_description();
    assert_eq!(description.reduced_degree(), Some(1));
    assert_eq!(description.infinitesimal_degree(), Some(41 * 41));
    assert_eq!(description.degree(), Some(41 * 41));
    assert_eq!(
        description.rational_points(),
        Some([isogeny.domain().identity()].as_slice())
    );
    assert!(matches!(description, KernelDescription::Mixed(_)));
}

#[test]
fn scalar_characteristic_factorization_splits_off_the_prime_to_p_part() {
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve(), 82).expect("scalar isogeny should build");
    let factorization = isogeny.scalar_characteristic_factorization();

    assert_eq!(factorization.p_power_exponent(), 1);
    assert_eq!(
        factorization.separable_part(),
        &num_bigint::BigUint::from(2u8)
    );
    assert_eq!(factorization.separable_degree(), 4);
    assert_eq!(factorization.infinitesimal_degree(), 41 * 41);
}

#[test]
fn visible_reduced_kernel_points_for_characteristic_divisible_scalar_come_from_the_prime_to_p_part()
{
    let curve = curve();
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 82).expect("scalar isogeny should build");

    let mut expected = vec![curve.identity()];
    expected.extend(curve.points_of_order(2));

    assert_eq!(
        isogeny
            .visible_reduced_kernel_points()
            .expect("visible reduced points should enumerate"),
        expected
    );
}

#[test]
fn mixed_kernel_description_for_p_times_m_uses_the_visible_m_torsion_and_residual_p_power_degree() {
    let curve = curve();
    let isogeny =
        ScalarMultiplicationIsogeny::new(curve.clone(), 82).expect("scalar isogeny should build");
    let description = isogeny.kernel_description();

    let mut expected = vec![curve.identity()];
    expected.extend(curve.points_of_order(2));

    assert_eq!(description.reduced_degree(), Some(4));
    assert_eq!(description.infinitesimal_degree(), Some(41 * 41));
    assert_eq!(description.degree(), Some(4 * 41 * 41));
    assert_eq!(description.rational_points(), Some(expected.as_slice()));
    match description {
        KernelDescription::Mixed(mixed) => {
            assert_eq!(
                mixed.label(),
                Some("kernel contribution from [n] = [p^1] o [2]")
            );
        }
        other => panic!("expected mixed kernel description, got {other:?}"),
    }
}
