use proptest::prelude::*;

use crate::elliptic_curves::{EndomorphismRingReport, ShortWeierstrassCurve};
use crate::fields::Fp;
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::frobenius::arb_frobenius_curve_case;

/// Endomorphism-side case derived from one finite short-Weierstrass curve.
#[derive(Clone, Debug)]
pub struct EndomorphismReportCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub report: EndomorphismRingReport,
}

/// Returns a finite-field curve together with its current endomorphism-ring
/// report.
pub fn arb_endomorphism_report_case<const P: u64>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<EndomorphismReportCase<P>> {
    arb_frobenius_curve_case::<P>(config)
        .prop_map(|case| EndomorphismReportCase {
            curve: case.curve,
            report: case
                .discriminant
                .endomorphism_ring_report()
                .expect("finite Frobenius discriminants should build a report"),
        })
        .boxed()
}

pub(crate) fn touch_endomorphism_case_fields() {
    let _ = |case: EndomorphismReportCase<17>| {
        let _ = case.curve;
        let _ = case.report;
    };
}
