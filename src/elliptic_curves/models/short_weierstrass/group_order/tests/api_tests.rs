use super::shared::{F43, f241_curve};
use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::group_order::{
        FiniteFieldGroupOrderStrategy, MestreConfig, SmallFieldSampledGroupOrderStrategy,
    },
};
use crate::fields::extension_field::define_fp_quadratic_extension;
use crate::fields::traits::Field;

define_fp_quadratic_extension!(
    spec: F43Sqrt2MestreSpec,
    field: F43Sqrt2Mestre,
    base: F43,
    non_residue: 2,
    name: "F43(sqrt(2)) for Mestre api tests",
);

#[test]
fn mestre_route_rejects_prime_fields_below_the_theorem_threshold() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let mut sampler = |_upper_bound: usize| Some(0usize);

    assert_eq!(
        curve.group_order_by_small_field_with_sampler(
            SmallFieldSampledGroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
            &mut sampler,
        ),
        Err(CurveError::MestrePrimeTooSmall { characteristic: 43 })
    );
}

#[test]
fn mestre_route_rejects_extension_fields() {
    let curve =
        ShortWeierstrassCurve::<F43Sqrt2Mestre>::new(F43Sqrt2Mestre::one(), F43Sqrt2Mestre::one())
            .expect("valid extension-field curve");
    let mut sampler = |_upper_bound: usize| Some(0usize);

    assert_eq!(
        curve.group_order_by_small_field_with_sampler(
            SmallFieldSampledGroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
            &mut sampler,
        ),
        Err(CurveError::MestreRequiresPrimeField {
            extension_degree: 2
        })
    );
}

#[test]
fn mestre_route_reports_iteration_cap_reached_before_sampling() {
    let curve = f241_curve();
    let mut sampler = |_upper_bound: usize| Some(0usize);

    assert_eq!(
        curve.group_order_by_small_field_with_sampler(
            SmallFieldSampledGroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(0)),
            &mut sampler,
        ),
        Err(CurveError::MestreIterationCapReached { max_iterations: 0 })
    );
}

#[test]
fn mestre_route_reports_sampler_exhaustion() {
    let curve = f241_curve();
    let mut sampler = |_upper_bound: usize| None::<usize>;

    assert_eq!(
        curve.group_order_by_small_field_with_sampler(
            SmallFieldSampledGroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
            &mut sampler,
        ),
        Err(CurveError::MestreSamplerExhausted)
    );
}

#[test]
fn finite_field_auto_route_uses_schoof_semantics() {
    let curve = ShortWeierstrassCurve::<crate::fields::Fp<7>>::new(
        crate::fields::Fp::<7>::from_i64(2),
        crate::fields::Fp::<7>::from_i64(3),
    )
    .expect("valid curve");

    let auto = curve
        .group_order_by(FiniteFieldGroupOrderStrategy::Auto)
        .expect("finite-field auto route should succeed");
    let schoof = curve
        .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
        .expect("explicit finite-field Schoof route should match finite-field auto");

    assert_eq!(auto, schoof);
}
