use crate::fields::traits::*;
use proptest::prelude::*;
use std::fmt;

use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::endomorphisms::EndomorphismRingReport;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::{CurveModel, FrobeniusTraceCurveModel};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::frobenius::arb_frobenius_curve_case;

/// Endomorphism-side case derived from one finite short-Weierstrass curve.
#[derive(Clone)]
pub struct EndomorphismReportCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub report: EndomorphismRingReport,
}

impl<F: Field> fmt::Debug for EndomorphismReportCase<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_struct("EndomorphismReportCase").finish()
    }
}

/// Returns a finite-field curve together with its current endomorphism-ring
/// report.
pub fn arb_endomorphism_report_case<F>(
    config: CurveStrategyConfig,
) -> BoxedStrategy<EndomorphismReportCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: 'static,
    ShortWeierstrassCurve<F>: CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>>
        + FrobeniusTraceCurveModel
        + 'static,
{
    arb_frobenius_curve_case::<F>(config)
        .prop_map(|case| {
            let report = EndomorphismRingReport::from_frobenius_trace(case.trace.clone())
                .expect("finite Frobenius discriminants should build a report");

            EndomorphismReportCase {
                curve: case.curve,
                report,
            }
        })
        .boxed()
}

pub(crate) fn touch_endomorphism_case_fields() {
    let _ = |case: EndomorphismReportCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.report;
    };
}
