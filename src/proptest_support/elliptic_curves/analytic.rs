use proptest::prelude::*;

use crate::elliptic_curves::{
    AnalyticWeierstrassCurve, ApproxTolerance, ComplexLattice, FundamentalParallelogramCoordinate,
    UpperHalfPlanePoint, WeierstrassCubicRoots,
};
use crate::proptest_support::config::AnalyticStrategyConfig;

/// Returns a validated upper-half-plane parameter `τ`.
pub fn arb_upper_half_plane_point(
    config: AnalyticStrategyConfig,
) -> BoxedStrategy<UpperHalfPlanePoint> {
    (
        -config.max_real_part..=config.max_real_part,
        config.min_imaginary_part..=config.max_imaginary_part,
    )
        .prop_map(|(re, im)| {
            UpperHalfPlanePoint::from_re_im(re, im)
                .expect("positive imaginary parts should stay in the upper half-plane")
        })
        .boxed()
}

/// Returns the standard lattice `Z + Zτ` attached to a sampled upper-half-plane
/// parameter.
pub fn arb_complex_lattice(config: AnalyticStrategyConfig) -> BoxedStrategy<ComplexLattice> {
    arb_upper_half_plane_point(config)
        .prop_map(ComplexLattice::from_tau)
        .boxed()
}

/// Returns a coordinate in the half-open unit square `[0,1) x [0,1)`.
pub fn arb_fundamental_coordinate() -> BoxedStrategy<FundamentalParallelogramCoordinate> {
    (0.0f64..1.0, 0.0f64..1.0)
        .prop_map(|(u, v)| {
            FundamentalParallelogramCoordinate::new(u, v)
                .expect("strategy stays inside the half-open unit square")
        })
        .boxed()
}

/// Returns a coordinate safely inside the unit square, away from the boundary.
pub fn arb_interior_fundamental_coordinate() -> BoxedStrategy<FundamentalParallelogramCoordinate> {
    (0.15f64..0.85, 0.15f64..0.85)
        .prop_map(|(u, v)| {
            FundamentalParallelogramCoordinate::new(u, v)
                .expect("interior strategy stays inside the half-open unit square")
        })
        .boxed()
}

/// Returns a stable non-singular analytic short-Weierstrass curve whose cubic
/// has three well-separated real roots summing to zero.
pub fn arb_stable_real_split_analytic_curve() -> BoxedStrategy<AnalyticWeierstrassCurve> {
    (0.4f64..3.0, 0.4f64..3.0)
        .prop_filter("real roots should stay well separated", |(e1, e2)| {
            let e3 = -(*e1 + *e2);
            (e1 - e2).abs() >= 0.2 && (e1 - e3).abs() >= 0.2 && (e2 - e3).abs() >= 0.2
        })
        .prop_map(|(e1, e2)| {
            let roots = WeierstrassCubicRoots::new(
                num_complex::Complex64::new(e1, 0.0),
                num_complex::Complex64::new(e2, 0.0),
                num_complex::Complex64::new(-(e1 + e2), 0.0),
                ApproxTolerance::strict(),
            )
            .expect("strategy only yields distinct real roots");
            AnalyticWeierstrassCurve::new(roots.g2(), roots.g3())
                .expect("distinct real roots should define a nonsingular curve")
        })
        .boxed()
}
