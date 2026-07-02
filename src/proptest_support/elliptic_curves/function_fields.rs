use crate::fields::traits::*;
use core::fmt;

use proptest::prelude::*;

use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::function_fields::{
    ShortWeierstrassFunction, ShortWeierstrassFunctionField,
};
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::proptest_support::config::{CurveStrategyConfig, PolynomialStrategyConfig};
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;
use crate::proptest_support::fields::arb_rational_function;

/// Ambient short-Weierstrass function field together with one sampled element.
pub struct FunctionFieldCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub field: ShortWeierstrassFunctionField<F>,
    pub function: ShortWeierstrassFunction<F>,
}

/// Ambient short-Weierstrass function field together with two sampled elements
/// on the same curve.
pub struct FunctionFieldPairCase<F: Field> {
    pub curve: ShortWeierstrassCurve<F>,
    pub field: ShortWeierstrassFunctionField<F>,
    pub left: ShortWeierstrassFunction<F>,
    pub right: ShortWeierstrassFunction<F>,
}

impl<F: Field> fmt::Debug for FunctionFieldCase<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FunctionFieldCase")
            .field("curve_a", self.curve.a())
            .field("curve_b", self.curve.b())
            .field("function", &self.function)
            .finish()
    }
}

impl<F: Field> fmt::Debug for FunctionFieldPairCase<F> {
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
pub fn arb_short_weierstrass_function_case<F>(
    curve_config: CurveStrategyConfig,
    polynomial_config: PolynomialStrategyConfig,
) -> BoxedStrategy<FunctionFieldCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(move |curve| {
            let a_strategy = arb_rational_function::<F>(polynomial_config);
            let b_strategy = arb_rational_function::<F>(polynomial_config);

            (Just(curve), a_strategy, b_strategy)
        })
        .prop_map(|(curve, a_part, b_part)| {
            let field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(curve.clone());
            let function = ShortWeierstrassFunction::<F>::new(curve.clone(), a_part, b_part);

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
pub fn arb_short_weierstrass_function_pair_case<F>(
    curve_config: CurveStrategyConfig,
    polynomial_config: PolynomialStrategyConfig,
) -> BoxedStrategy<FunctionFieldPairCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(move |curve| {
            let left_a = arb_rational_function::<F>(polynomial_config);
            let left_b = arb_rational_function::<F>(polynomial_config);
            let right_a = arb_rational_function::<F>(polynomial_config);
            let right_b = arb_rational_function::<F>(polynomial_config);

            (Just(curve), left_a, left_b, right_a, right_b)
        })
        .prop_map(|(curve, left_a, left_b, right_a, right_b)| {
            let field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(curve.clone());
            let left = ShortWeierstrassFunction::<F>::new(curve.clone(), left_a, left_b);
            let right = ShortWeierstrassFunction::<F>::new(curve.clone(), right_a, right_b);

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
    let _ = |case: FunctionFieldCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.field;
        let _ = case.function;
    };
    let _ = |case: FunctionFieldPairCase<crate::fields::Fp17>| {
        let _ = case.curve;
        let _ = case.field;
        let _ = case.left;
        let _ = case.right;
    };
}
