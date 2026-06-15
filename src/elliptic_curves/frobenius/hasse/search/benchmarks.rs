use std::hint::black_box;
use std::time::Instant;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    frobenius::{
        HasseInterval,
        hasse::search::{HasseBsgsConfig, HasseBsgsParity, HasseBsgsTraversal},
    },
    models::short_weierstrass::group_order_parity::GroupOrderParity,
    traits::{AffineCurveModel, HasseIntervalSearchCurveModel},
};
use crate::fields::{Fp, traits::Field};

type FLarge = Fp<1_000_000_007>;

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing even parity against unknown parity"]
fn benchmark_even_parity_vs_unknown_bsgs() {
    let curve =
        ShortWeierstrassCurve::<FLarge>::new(FLarge::zero(), FLarge::one()).expect("valid curve");
    let point = curve
        .point(FLarge::zero(), FLarge::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval =
        HasseInterval::for_q(1_000_000_007).expect("benchmark Hasse interval should exist");
    let repetitions = 40_000usize;

    assert_eq!(
        curve.group_order_parity_from_two_torsion(),
        GroupOrderParity::Even
    );

    let unknown_config = HasseBsgsConfig::new();
    let even_config = HasseBsgsConfig::new().with_known_parity(HasseBsgsParity::Even);

    for _ in 0..64 {
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                unknown_config,
            ),
        );
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                even_config,
            ),
        );
    }

    let unknown_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, unknown_config, repetitions);
    let even_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, even_config, repetitions);

    eprintln!(
        "unknown={:?} even={:?} reps={repetitions} unknown_ns_per_iter={} even_ns_per_iter={}",
        unknown_elapsed,
        even_elapsed,
        unknown_elapsed.as_nanos() / repetitions as u128,
        even_elapsed.as_nanos() / repetitions as u128
    );
}

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing odd parity against unknown parity"]
fn benchmark_odd_parity_vs_unknown_bsgs() {
    let curve = ShortWeierstrassCurve::<FLarge>::new(FLarge::from_i64(3), FLarge::one())
        .expect("valid curve");
    let point = curve
        .point(FLarge::zero(), FLarge::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval = crate::elliptic_curves::frobenius::HasseInterval::for_q(1_000_000_007)
        .expect("benchmark Hasse interval should exist");
    let repetitions = 40_000usize;

    assert_eq!(
        curve.group_order_parity_from_two_torsion(),
        GroupOrderParity::Odd
    );

    let unknown_config = HasseBsgsConfig::new();
    let odd_config = HasseBsgsConfig::new().with_known_parity(HasseBsgsParity::Odd);

    for _ in 0..64 {
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                unknown_config,
            ),
        );
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                odd_config,
            ),
        );
    }

    let unknown_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, unknown_config, repetitions);
    let odd_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, odd_config, repetitions);

    eprintln!(
        "unknown={:?} odd={:?} reps={repetitions} unknown_ns_per_iter={} odd_ns_per_iter={}",
        unknown_elapsed,
        odd_elapsed,
        unknown_elapsed.as_nanos() / repetitions as u128,
        odd_elapsed.as_nanos() / repetitions as u128
    );
}

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing BSGS traversals"]
fn benchmark_middle_out_vs_left_to_right() {
    // This is the fixed-instance benchmark. It answers:
    // "on one concrete curve/point pair, does changing only the traversal
    // order help or hurt?"
    //
    // It is intentionally *not* a benchmark for the semicircular heuristic as
    // a distributional claim; for that we use the center-heavy corpus below.
    let curve = ShortWeierstrassCurve::<FLarge>::new(FLarge::one(), FLarge::one())
        .expect("valid benchmark curve");
    let point = curve
        .point(FLarge::zero(), FLarge::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval =
        HasseInterval::for_q(1_000_000_007).expect("benchmark Hasse interval should exist");
    let repetitions = 40_000usize;

    let left_to_right_config = HasseBsgsConfig::new()
        .with_traversal(HasseBsgsTraversal::LeftToRight)
        .with_fast_negation(true);
    let middle_out_config = HasseBsgsConfig::new()
        .with_traversal(HasseBsgsTraversal::MiddleOut)
        .with_fast_negation(true);

    for _ in 0..64 {
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                left_to_right_config,
            ),
        );
        let _ = black_box(
            curve.find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                middle_out_config,
            ),
        );
    }

    let left_to_right_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, left_to_right_config, repetitions);
    let middle_out_elapsed =
        benchmark_bsgs_configuration(&curve, &point, &interval, middle_out_config, repetitions);

    eprintln!(
        "left_to_right={:?} middle_out={:?} reps={repetitions} left_ns_per_iter={} middle_ns_per_iter={}",
        left_to_right_elapsed,
        middle_out_elapsed,
        left_to_right_elapsed.as_nanos() / repetitions as u128,
        middle_out_elapsed.as_nanos() / repetitions as u128
    );
}

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing traversals on a center-heavy corpus"]
fn benchmark_middle_out_vs_left_to_right_on_center_heavy_corpus() {
    // This benchmark is designed to let the center-first heuristic show its
    // intended behavior. Instead of one fixed instance, it uses a small
    // deterministic corpus of instances whose decisive annihilating multiple
    // lies unusually close to the center `q + 1` of the Hasse interval.
    //
    // In other words:
    // - the fixed-instance benchmark measures raw traversal overhead
    // - this corpus benchmark measures the distributional regime that should
    //   favor MiddleOut
    let interval =
        HasseInterval::for_q(1_000_000_007).expect("benchmark Hasse interval should exist");
    let corpus = build_center_heavy_benchmark_corpus(&interval);
    let repetitions = 400usize;

    let left_to_right_config = HasseBsgsConfig::new()
        .with_traversal(HasseBsgsTraversal::LeftToRight)
        .with_fast_negation(true);
    let middle_out_config = HasseBsgsConfig::new()
        .with_traversal(HasseBsgsTraversal::MiddleOut)
        .with_fast_negation(true);

    for _ in 0..16 {
        let _ = black_box(benchmark_bsgs_corpus_iteration(
            &corpus,
            &interval,
            left_to_right_config,
        ));
        let _ = black_box(benchmark_bsgs_corpus_iteration(
            &corpus,
            &interval,
            middle_out_config,
        ));
    }

    let left_to_right_start = Instant::now();
    for _ in 0..repetitions {
        let checksum = benchmark_bsgs_corpus_iteration(&corpus, &interval, left_to_right_config);
        black_box(checksum);
    }
    let left_to_right_elapsed = left_to_right_start.elapsed();

    let middle_out_start = Instant::now();
    for _ in 0..repetitions {
        let checksum = benchmark_bsgs_corpus_iteration(&corpus, &interval, middle_out_config);
        black_box(checksum);
    }
    let middle_out_elapsed = middle_out_start.elapsed();

    let average_distance_to_center = corpus
        .iter()
        .map(|instance| distance_to_interval_center(instance.annihilating_multiple, &interval))
        .sum::<u128>()
        / corpus.len() as u128;

    eprintln!(
        "center_heavy_cases={} avg_center_distance={} left_to_right={:?} middle_out={:?} reps={repetitions} left_ns_per_iter={} middle_ns_per_iter={}",
        corpus.len(),
        average_distance_to_center,
        left_to_right_elapsed,
        middle_out_elapsed,
        left_to_right_elapsed.as_nanos() / repetitions as u128,
        middle_out_elapsed.as_nanos() / repetitions as u128
    );
}

fn benchmark_bsgs_configuration(
    curve: &ShortWeierstrassCurve<FLarge>,
    point: &crate::elliptic_curves::AffinePoint<FLarge>,
    interval: &crate::elliptic_curves::frobenius::HasseInterval,
    config: HasseBsgsConfig,
    repetitions: usize,
) -> std::time::Duration {
    let start = Instant::now();
    for _ in 0..repetitions {
        let result = curve
            .find_annihilating_multiple_in_interval_bsgs_with_config(
                point,
                interval.clone(),
                config,
            )
            .expect("BSGS benchmark configuration should succeed");
        black_box(result);
    }
    start.elapsed()
}

/// One deterministic benchmark instance for the center-heavy corpus.
///
/// The stored `annihilating_multiple` is not used by the runtime benchmark
/// loop itself; it is recorded so that corpus construction can sort instances
/// by distance to the center of `H(q)` and so that the printed summary can say
/// how “center-heavy” the final corpus really is.
#[derive(Clone)]
struct CenterHeavyBenchmarkInstance {
    curve: ShortWeierstrassCurve<FLarge>,
    point: crate::elliptic_curves::AffinePoint<FLarge>,
    annihilating_multiple: u128,
}

fn build_center_heavy_benchmark_corpus(
    interval: &HasseInterval,
) -> Vec<CenterHeavyBenchmarkInstance> {
    // We generate a deterministic family of curves with the same easy point
    // `(0, 1)`, recover one annihilating multiple on each curve by the
    // left-to-right baseline, then keep the instances whose annihilating
    // multiple is closest to `q + 1`.
    //
    // This makes the corpus reproducible while still approximating the regime
    // in which a center-first traversal should terminate earlier.
    let selection_target = 32usize;
    let search_budget = 512u64;
    let mut candidates = Vec::new();

    let left_to_right_config = HasseBsgsConfig::new()
        .with_traversal(HasseBsgsTraversal::LeftToRight)
        .with_fast_negation(true);

    for a in 1..=search_budget {
        let Ok(curve) =
            ShortWeierstrassCurve::<FLarge>::new(FLarge::from_i64(a as i64), FLarge::one())
        else {
            continue;
        };
        let point = curve
            .point(FLarge::zero(), FLarge::one())
            .expect("(0, 1) should lie on the benchmark curve family");
        let Some(annihilating_multiple) = curve
            .find_annihilating_multiple_in_interval_bsgs_with_config(
                &point,
                interval.clone(),
                left_to_right_config,
            )
            .expect("benchmark corpus search should succeed")
        else {
            continue;
        };

        candidates.push(CenterHeavyBenchmarkInstance {
            curve,
            point,
            annihilating_multiple,
        });
    }

    candidates.sort_by_key(|instance| {
        distance_to_interval_center(instance.annihilating_multiple, interval)
    });
    candidates.truncate(selection_target);
    assert!(
        !candidates.is_empty(),
        "center-heavy benchmark corpus should contain at least one instance"
    );
    candidates
}

/// Runs one full benchmark pass over the deterministic corpus.
///
/// The XOR checksum is only there to keep the optimizer honest without
/// polluting the benchmark with extra reporting work inside the hot loop.
fn benchmark_bsgs_corpus_iteration(
    corpus: &[CenterHeavyBenchmarkInstance],
    interval: &HasseInterval,
    config: HasseBsgsConfig,
) -> u128 {
    let mut checksum = 0u128;
    for instance in corpus {
        let annihilating_multiple = instance
            .curve
            .find_annihilating_multiple_in_interval_bsgs_with_config(
                &instance.point,
                interval.clone(),
                config,
            )
            .expect("benchmark corpus configuration should succeed")
            .expect("benchmark corpus search should find an annihilating multiple");
        checksum ^= annihilating_multiple;
    }
    checksum
}

/// Distance from one candidate to the center `q + 1` of the Hasse interval.
///
/// This is the geometric quantity that the middle-out heuristic cares about:
/// smaller values mean “closer to where the center-first traversal looks
/// first”.
fn distance_to_interval_center(candidate: u128, interval: &HasseInterval) -> u128 {
    let center = interval.q() + 1;
    candidate.abs_diff(center)
}
