use core::fmt;

use proptest::prelude::*;

use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::function_fields::{
    ShortWeierstrassFunction, ShortWeierstrassFunctionField,
};
use crate::fields::Fp;
use crate::proptest_support::config::{CurveStrategyConfig, PolynomialStrategyConfig};
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;
use crate::proptest_support::fields::arb_rational_function;

/// Ambient short-Weierstrass function field together with one sampled element.
pub struct FunctionFieldCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub field: ShortWeierstrassFunctionField<Fp<P>>,
    pub function: ShortWeierstrassFunction<Fp<P>>,
}

/// Ambient short-Weierstrass function field together with two sampled elements
/// on the same curve.
pub struct FunctionFieldPairCase<const P: u64> {
    pub curve: ShortWeierstrassCurve<Fp<P>>,
    pub field: ShortWeierstrassFunctionField<Fp<P>>,
    pub left: ShortWeierstrassFunction<Fp<P>>,
    pub right: ShortWeierstrassFunction<Fp<P>>,
}

impl<const P: u64> fmt::Debug for FunctionFieldCase<P> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FunctionFieldCase")
            .field("curve_a", self.curve.a())
            .field("curve_b", self.curve.b())
            .field("function", &self.function)
            .finish()
    }
}

impl<const P: u64> fmt::Debug for FunctionFieldPairCase<P> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FunctionFieldPairCase")
            .field("curve_a", self.curve.a())
            .field("curve_b", self.curve.b())
            .field("left", &self.left)
            .field("right", &self.right)
            .finish()
    }
}

/// Returns one nonsingular short-Weierstrass curve together with one sampled
/// function-field element over that same curve.
pub fn arb_short_weierstrass_function_case<const P: u64>(
    curve_config: CurveStrategyConfig,
    polynomial_config: PolynomialStrategyConfig,
) -> BoxedStrategy<FunctionFieldCase<P>> {
    arb_nonsingular_curve::<P>(curve_config)
        .prop_flat_map(move |curve| {
            let a_strategy = arb_rational_function::<Fp<P>>(polynomial_config);
            let b_strategy = arb_rational_function::<Fp<P>>(polynomial_config);

            (Just(curve), a_strategy, b_strategy)
        })
        .prop_map(|(curve, a_part, b_part)| {
            let field = ShortWeierstrassFunctionField::<Fp<P>>::new(curve.clone());
            let function = ShortWeierstrassFunction::<Fp<P>>::new(curve.clone(), a_part, b_part);

            FunctionFieldCase {
                curve,
                field,
                function,
            }
        })
        .boxed()
}

/// Returns one nonsingular short-Weierstrass curve together with two sampled
/// function-field elements over that same curve.
pub fn arb_short_weierstrass_function_pair_case<const P: u64>(
    curve_config: CurveStrategyConfig,
    polynomial_config: PolynomialStrategyConfig,
) -> BoxedStrategy<FunctionFieldPairCase<P>> {
    arb_nonsingular_curve::<P>(curve_config)
        .prop_flat_map(move |curve| {
            let left_a = arb_rational_function::<Fp<P>>(polynomial_config);
            let left_b = arb_rational_function::<Fp<P>>(polynomial_config);
            let right_a = arb_rational_function::<Fp<P>>(polynomial_config);
            let right_b = arb_rational_function::<Fp<P>>(polynomial_config);

            (Just(curve), left_a, left_b, right_a, right_b)
        })
        .prop_map(|(curve, left_a, left_b, right_a, right_b)| {
            let field = ShortWeierstrassFunctionField::<Fp<P>>::new(curve.clone());
            let left = ShortWeierstrassFunction::<Fp<P>>::new(curve.clone(), left_a, left_b);
            let right = ShortWeierstrassFunction::<Fp<P>>::new(curve.clone(), right_a, right_b);

            FunctionFieldPairCase {
                curve,
                field,
                left,
                right,
            }
        })
        .boxed()
}

pub(crate) fn touch_function_field_case_fields() {
    let _ = |case: FunctionFieldCase<17>| {
        let _ = case.curve;
        let _ = case.field;
        let _ = case.function;
    };
    let _ = |case: FunctionFieldPairCase<17>| {
        let _ = case.curve;
        let _ = case.field;
        let _ = case.left;
        let _ = case.right;
    };
}
