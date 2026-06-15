use std::hint::black_box;
use std::time::Instant;

use num_bigint::BigUint;

use super::shared::{bu, f7_curve, f7_point};

#[test]
#[ignore = "benchmark helper; run explicitly in release mode with --nocapture"]
fn benchmark_point_order_from_multiple_optimized_against_baseline() {
    let curve = f7_curve();
    let point = f7_point(6, 0);
    let multiple =
        BigUint::from(2u8).pow(20) * BigUint::from(3u8).pow(12) * BigUint::from(5u8).pow(8);
    let factorization = vec![(bu(2), 20), (bu(3), 12), (bu(5), 8)];
    let iterations = 2_000usize;

    let baseline_start = Instant::now();
    for _ in 0..iterations {
        let report = crate::elliptic_curves::models::short_weierstrass::point_order::point_order_from_multiple_baseline(
            &curve,
            &point,
            multiple.clone(),
            &factorization,
        )
        .expect("baseline report should build");
        black_box(report);
    }
    let baseline_elapsed = baseline_start.elapsed();

    let optimized_start = Instant::now();
    for _ in 0..iterations {
        let report = curve
            .point_order_from_multiple(&point, multiple.clone(), &factorization)
            .expect("optimized report should build");
        black_box(report);
    }
    let optimized_elapsed = optimized_start.elapsed();

    println!("iterations: {iterations}");
    println!("baseline:  {:?}", baseline_elapsed);
    println!("optimized: {:?}", optimized_elapsed);
    println!(
        "speedup:   {:.2}x",
        baseline_elapsed.as_secs_f64() / optimized_elapsed.as_secs_f64()
    );
}
