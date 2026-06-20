use std::hint::black_box;
use std::time::{Duration, Instant};

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::group_order::{
        FiniteFieldGroupOrderStrategy, MestreConfig, SmallFieldSampledGroupOrderStrategy,
    },
};
use crate::fields::{Fp, traits::Field};

use super::shared::{F241, genuine_twist_curve, max_order_point_index};

type FLarge = Fp<1_000_000_007>;

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing Schoof against Mestre"]
fn benchmark_schoof_vs_mestre_on_f241() {
    let (curve, original_index, twist_index) = find_schoof_and_mestre_benchmark_curve();
    let repetitions = 400usize;
    let mestre_config = MestreConfig::with_iteration_cap(8);

    for _ in 0..16 {
        let _ = black_box(
            curve
                .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
                .expect("benchmark curve should resolve under Schoof"),
        );
        let _ = black_box(run_mestre_once(
            &curve,
            mestre_config.clone(),
            original_index,
            twist_index,
        ));
    }

    let schoof_elapsed = benchmark_fixed_group_order_strategy(
        &curve,
        FiniteFieldGroupOrderStrategy::Schoof,
        repetitions,
    );
    let mestre_elapsed = benchmark_mestre_strategy(
        &curve,
        mestre_config,
        original_index,
        twist_index,
        repetitions,
    );

    eprintln!(
        "curve_a={:?} curve_b={:?} schoof={:?} mestre={:?} reps={repetitions} schoof_ns_per_iter={} mestre_ns_per_iter={}",
        curve.a(),
        curve.b(),
        schoof_elapsed,
        mestre_elapsed,
        schoof_elapsed.as_nanos() / repetitions as u128,
        mestre_elapsed.as_nanos() / repetitions as u128
    );
}

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when benchmarking Schoof over Fp<10^9 + 7>"]
fn benchmark_schoof_on_fp_1e9_plus_7() {
    let curve = ShortWeierstrassCurve::<FLarge>::new(FLarge::from_i64(-8), FLarge::zero())
        .expect("benchmark curve y^2 = x^3 - 8x over Fp<10^9 + 7> should be nonsingular");
    let repetitions = 20usize;

    for _ in 0..4 {
        let _ = black_box(
            curve
                .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
                .expect("benchmark curve should resolve under Schoof"),
        );
    }

    let schoof_elapsed = benchmark_fixed_group_order_strategy_large_prime(
        &curve,
        FiniteFieldGroupOrderStrategy::Schoof,
        repetitions,
    );

    eprintln!(
        "curve_a={:?} curve_b={:?} schoof={:?} reps={repetitions} schoof_ns_per_iter={}",
        curve.a(),
        curve.b(),
        schoof_elapsed,
        schoof_elapsed.as_nanos() / repetitions as u128,
    );
}

fn find_schoof_and_mestre_benchmark_curve() -> (ShortWeierstrassCurve<F241>, usize, usize) {
    for a in -8..=8 {
        for b in -8..=8 {
            let Ok(curve) =
                ShortWeierstrassCurve::<F241>::new(F241::from_i64(a), F241::from_i64(b))
            else {
                continue;
            };
            if curve
                .group_order_by(FiniteFieldGroupOrderStrategy::Schoof)
                .is_err()
            {
                continue;
            }

            let twist_curve = genuine_twist_curve(&curve);
            let original_index = max_order_point_index(&curve);
            let twist_index = max_order_point_index(&twist_curve);
            if run_mestre_once(
                &curve,
                MestreConfig::with_iteration_cap(8),
                original_index,
                twist_index,
            )
            .is_ok()
            {
                return (curve, original_index, twist_index);
            }
        }
    }

    panic!("expected to find one deterministic F241 curve resolved by both Schoof and Mestre");
}

fn benchmark_fixed_group_order_strategy(
    curve: &ShortWeierstrassCurve<F241>,
    strategy: FiniteFieldGroupOrderStrategy,
    repetitions: usize,
) -> Duration {
    let start = Instant::now();
    for _ in 0..repetitions {
        let report = curve
            .group_order_by(strategy)
            .expect("benchmark strategy should succeed");
        black_box(report);
    }
    start.elapsed()
}

fn benchmark_fixed_group_order_strategy_large_prime(
    curve: &ShortWeierstrassCurve<FLarge>,
    strategy: FiniteFieldGroupOrderStrategy,
    repetitions: usize,
) -> Duration {
    let start = Instant::now();
    for _ in 0..repetitions {
        let report = curve
            .group_order_by(strategy)
            .expect("benchmark strategy should succeed");
        black_box(report);
    }
    start.elapsed()
}

fn benchmark_mestre_strategy(
    curve: &ShortWeierstrassCurve<F241>,
    config: MestreConfig,
    original_index: usize,
    twist_index: usize,
    repetitions: usize,
) -> Duration {
    let start = Instant::now();
    for _ in 0..repetitions {
        let report = run_mestre_once(curve, config.clone(), original_index, twist_index)
            .expect("benchmark Mestre run should succeed");
        black_box(report);
    }
    start.elapsed()
}

fn run_mestre_once(
    curve: &ShortWeierstrassCurve<F241>,
    config: MestreConfig,
    original_index: usize,
    twist_index: usize,
) -> Result<
    crate::elliptic_curves::frobenius::group_order::GroupOrderReport,
    crate::elliptic_curves::CurveError,
> {
    let mut requested = vec![original_index, twist_index].into_iter();
    let mut sampler = move |_upper_bound: usize| requested.next().or(Some(original_index));
    curve.group_order_by_small_field_with_sampler(
        SmallFieldSampledGroupOrderStrategy::MestreFp(config),
        &mut sampler,
    )
}
