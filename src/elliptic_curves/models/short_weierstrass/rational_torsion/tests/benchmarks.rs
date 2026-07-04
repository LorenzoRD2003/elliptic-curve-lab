use std::hint::black_box;
use std::time::{Duration, Instant};

use proptest::prelude::*;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::rational_torsion::{RationalTorsionReport, RationalTorsionStrategy},
};
use crate::fields::Q;

use super::fixtures::{
    cyclic_six_fixture, product_two_two_fixture, q, rational_scaled_fixture,
    trivial_torsion_fixture,
};

const BENCHMARK_REPETITIONS: usize = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FasterStrategy {
    LutzNagell,
    GoodReductionHensel,
    Tie,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct StrategyTiming {
    strategy: RationalTorsionStrategy,
    elapsed: Duration,
    report_count: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct StrategyBenchmarkConclusion {
    lutz_nagell: StrategyTiming,
    good_reduction_hensel: StrategyTiming,
    faster_strategy: FasterStrategy,
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn rational_torsion_strategies_agree_on_small_integral_curves(
        a in -4i64..=4,
        b in -4i64..=4,
    ) {
        let Ok(curve) = integral_curve(a, b) else {
            prop_assume!(false);
            return Ok(());
        };

        let lutz_nagell = report_by(&curve, RationalTorsionStrategy::LutzNagell);
        let good_reduction_hensel = report_by(&curve, RationalTorsionStrategy::GoodReductionHensel);

        prop_assert!(
            lutz_nagell.is_ok(),
            "Lutz-Nagell should classify small integral fixture y² = x³ + {a}x + {b}: {lutz_nagell:?}"
        );
        prop_assert!(
            good_reduction_hensel.is_ok(),
            "good-reduction/Hensel should classify small integral fixture y² = x³ + {a}x + {b}: {good_reduction_hensel:?}"
        );

        assert_same_torsion_payload(&lutz_nagell.unwrap(), &good_reduction_hensel.unwrap());
    }
}

#[test]
#[ignore = "timing benchmark; run with --ignored --nocapture when comparing strategy speed"]
fn benchmark_rational_torsion_strategies_on_small_curve_corpus() {
    let corpus = benchmark_corpus();
    assert!(
        !corpus.is_empty(),
        "benchmark corpus should contain at least one non-singular curve"
    );
    assert_strategies_agree_on_corpus(&corpus);

    let lutz_nagell = time_strategy(&corpus, RationalTorsionStrategy::LutzNagell);
    let good_reduction_hensel =
        time_strategy(&corpus, RationalTorsionStrategy::GoodReductionHensel);
    let conclusion = StrategyBenchmarkConclusion::new(lutz_nagell, good_reduction_hensel);

    println!("{}", conclusion.describe());
}

impl StrategyBenchmarkConclusion {
    fn new(lutz_nagell: StrategyTiming, good_reduction_hensel: StrategyTiming) -> Self {
        let faster_strategy = match lutz_nagell.elapsed.cmp(&good_reduction_hensel.elapsed) {
            std::cmp::Ordering::Less => FasterStrategy::LutzNagell,
            std::cmp::Ordering::Equal => FasterStrategy::Tie,
            std::cmp::Ordering::Greater => FasterStrategy::GoodReductionHensel,
        };

        Self {
            lutz_nagell,
            good_reduction_hensel,
            faster_strategy,
        }
    }

    fn describe(&self) -> String {
        let speedup = speedup_ratio(self.lutz_nagell.elapsed, self.good_reduction_hensel.elapsed);
        format!(
            "\
rational torsion strategy benchmark
corpus reports per strategy: {}
Lutz-Nagell: {:?}
good-reduction/Hensel: {:?}
faster strategy: {:?}
speedup ratio: {:.2}×",
            self.lutz_nagell.report_count,
            self.lutz_nagell.elapsed,
            self.good_reduction_hensel.elapsed,
            self.faster_strategy,
            speedup
        )
    }
}

fn benchmark_corpus() -> Vec<ShortWeierstrassCurve<Q>> {
    [
        product_two_two_fixture(),
        cyclic_six_fixture(),
        trivial_torsion_fixture(),
        rational_scaled_fixture(),
    ]
    .into_iter()
    .map(|fixture| fixture.curve)
    .collect()
}

fn integral_curve(
    a: i64,
    b: i64,
) -> Result<ShortWeierstrassCurve<Q>, crate::elliptic_curves::CurveError> {
    ShortWeierstrassCurve::<Q>::new(q(a, 1), q(b, 1))
}

fn assert_strategies_agree_on_corpus(corpus: &[ShortWeierstrassCurve<Q>]) {
    for curve in corpus {
        let lutz_nagell = report_by(curve, RationalTorsionStrategy::LutzNagell)
            .expect("Lutz-Nagell should classify benchmark curve");
        let good_reduction_hensel = report_by(curve, RationalTorsionStrategy::GoodReductionHensel)
            .expect("good-reduction/Hensel should classify benchmark curve");
        assert_same_torsion_payload(&lutz_nagell, &good_reduction_hensel);
    }
}

fn assert_same_torsion_payload(left: &RationalTorsionReport, right: &RationalTorsionReport) {
    assert_eq!(left.group(), right.group());
    assert_eq!(left.points(), right.points());
}

fn time_strategy(
    corpus: &[ShortWeierstrassCurve<Q>],
    strategy: RationalTorsionStrategy,
) -> StrategyTiming {
    let started_at = Instant::now();
    let mut report_count = 0usize;

    for _ in 0..BENCHMARK_REPETITIONS {
        for curve in corpus {
            let report = report_by(black_box(curve), strategy)
                .expect("benchmark corpus should be classifiable by both strategies");
            black_box(report);
            report_count += 1;
        }
    }

    StrategyTiming {
        strategy,
        elapsed: started_at.elapsed(),
        report_count,
    }
}

fn report_by(
    curve: &ShortWeierstrassCurve<Q>,
    strategy: RationalTorsionStrategy,
) -> Result<
    RationalTorsionReport,
    crate::elliptic_curves::short_weierstrass::rational_torsion::RationalTorsionError,
> {
    curve.rational_torsion_by(strategy)
}

fn speedup_ratio(lutz_nagell: Duration, good_reduction_hensel: Duration) -> f64 {
    let lutz = lutz_nagell.as_secs_f64();
    let hensel = good_reduction_hensel.as_secs_f64();
    if lutz == 0.0 || hensel == 0.0 {
        return 1.0;
    }
    lutz.max(hensel) / lutz.min(hensel)
}
