use crate::elliptic_curves::frobenius::hasse::search_bsgs::{
    HasseBsgsConfig, HasseBsgsParity, HasseBsgsTraversal,
    find_annihilating_multiple_in_interval_bsgs_with_config,
};
use crate::elliptic_curves::{AffineCurveModel, ShortWeierstrassCurve};
use crate::fields::{Field, Fp};

use std::hint::black_box;
use std::time::Instant;

type FLarge = Fp<1_000_000_007>;

#[test]
#[ignore = "microbenchmark; run explicitly in release mode when comparing BSGS configurations"]
fn benchmark_fast_negation_vs_plain_bsgs() {
    let curve = ShortWeierstrassCurve::<FLarge>::new(FLarge::one(), FLarge::one())
        .expect("valid benchmark curve");
    let point = curve
        .point(FLarge::zero(), FLarge::one())
        .expect("(0, 1) should lie on the benchmark curve");
    let interval = crate::elliptic_curves::HasseInterval::for_q(1_000_000_007)
        .expect("benchmark Hasse interval should exist");
    let repetitions = 40_000usize;

    let plain_config = HasseBsgsConfig {
        traversal: HasseBsgsTraversal::LeftToRight,
        use_fast_negation: false,
        known_parity: HasseBsgsParity::Unknown,
    };
    let fast_config = HasseBsgsConfig {
        traversal: HasseBsgsTraversal::LeftToRight,
        use_fast_negation: true,
        known_parity: HasseBsgsParity::Unknown,
    };

    for _ in 0..256 {
        let _ = black_box(find_annihilating_multiple_in_interval_bsgs_with_config(
            &curve,
            &point,
            interval.clone(),
            plain_config,
        ));
        let _ = black_box(find_annihilating_multiple_in_interval_bsgs_with_config(
            &curve,
            &point,
            interval.clone(),
            fast_config,
        ));
    }

    let plain_start = Instant::now();
    for _ in 0..repetitions {
        let result = find_annihilating_multiple_in_interval_bsgs_with_config(
            &curve,
            &point,
            interval.clone(),
            plain_config,
        )
        .expect("plain BSGS should succeed");
        black_box(result);
    }
    let plain_elapsed = plain_start.elapsed();

    let fast_start = Instant::now();
    for _ in 0..repetitions {
        let result = find_annihilating_multiple_in_interval_bsgs_with_config(
            &curve,
            &point,
            interval.clone(),
            fast_config,
        )
        .expect("fast-negation BSGS should succeed");
        black_box(result);
    }
    let fast_elapsed = fast_start.elapsed();

    eprintln!(
        "plain={:?} fast={:?} reps={repetitions} plain_ns_per_iter={} fast_ns_per_iter={}",
        plain_elapsed,
        fast_elapsed,
        plain_elapsed.as_nanos() / repetitions as u128,
        fast_elapsed.as_nanos() / repetitions as u128
    );
}
