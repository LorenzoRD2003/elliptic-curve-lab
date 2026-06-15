use super::shared::{F43, f241_curve};
use crate::elliptic_curves::{
    CurveError, ShortWeierstrassCurve,
    frobenius::group_order::{GroupOrderStrategy, MestreConfig},
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
fn deterministic_group_order_api_reports_that_mestre_needs_a_sampler() {
    let curve = f241_curve();

    assert_eq!(
        curve.group_order_by(GroupOrderStrategy::MestreFp(MestreConfig::unbounded())),
        Err(CurveError::GroupOrderStrategyRequiresSampler {
            strategy: "MestreFp"
        })
    );
}

#[test]
fn mestre_route_rejects_prime_fields_below_the_theorem_threshold() {
    let curve = ShortWeierstrassCurve::<F43>::new(F43::one(), F43::one()).expect("valid curve");
    let mut sampler = |_upper_bound: usize| Some(0usize);

    assert_eq!(
        curve.group_order_by_with_sampler(
            GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
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
        curve.group_order_by_with_sampler(
            GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
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
        curve.group_order_by_with_sampler(
            GroupOrderStrategy::MestreFp(MestreConfig::with_iteration_cap(0)),
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
        curve.group_order_by_with_sampler(
            GroupOrderStrategy::MestreFp(MestreConfig::unbounded()),
            &mut sampler,
        ),
        Err(CurveError::MestreSamplerExhausted)
    );
}
