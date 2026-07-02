use crate::fields::traits::*;
use core::fmt;

use proptest::prelude::*;

use crate::elliptic_curves::AffinePoint;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::short_weierstrass::function_fields::{
    ShortWeierstrassFunction, ShortWeierstrassFunctionField,
};
use crate::elliptic_curves::short_weierstrass::isogenies::function_field_maps::ShortWeierstrassFunctionFieldMap;
use crate::elliptic_curves::traits::{CurveModel, EnumerableCurveModel};
use crate::fields::{
    rational_function_field::RationalFunction,
    traits::{EnumerableFiniteField, SqrtField},
};
use crate::proptest_support::config::CurveStrategyConfig;
use crate::proptest_support::elliptic_curves::short_weierstrass::arb_nonsingular_curve;

/// One valid short-Weierstrass function-field pullback map together with its
/// ambient domain and codomain function fields.
pub struct FunctionFieldMapCase<F: Field> {
    pub domain_curve: ShortWeierstrassCurve<F>,
    pub codomain_curve: ShortWeierstrassCurve<F>,
    pub domain_field: ShortWeierstrassFunctionField<F>,
    pub codomain_field: ShortWeierstrassFunctionField<F>,
    pub map: ShortWeierstrassFunctionFieldMap<F>,
}

/// Two composable short-Weierstrass function-field pullback maps together with
/// their computed composition.
pub struct ComposableFunctionFieldMapCase<F: Field> {
    pub domain_curve: ShortWeierstrassCurve<F>,
    pub middle_curve: ShortWeierstrassCurve<F>,
    pub codomain_curve: ShortWeierstrassCurve<F>,
    pub first: ShortWeierstrassFunctionFieldMap<F>,
    pub second: ShortWeierstrassFunctionFieldMap<F>,
    pub composite: ShortWeierstrassFunctionFieldMap<F>,
}

impl<F> fmt::Debug for FunctionFieldMapCase<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FunctionFieldMapCase")
            .field("domain_curve_a", self.domain_curve.a())
            .field("domain_curve_b", self.domain_curve.b())
            .field("codomain_curve_a", self.codomain_curve.a())
            .field("codomain_curve_b", self.codomain_curve.b())
            .field("x_pullback", self.map.x_pullback())
            .field("y_pullback", self.map.y_pullback())
            .finish()
    }
}

impl<F> fmt::Debug for ComposableFunctionFieldMapCase<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ComposableFunctionFieldMapCase")
            .field("domain_curve_a", self.domain_curve.a())
            .field("middle_curve_a", self.middle_curve.a())
            .field("codomain_curve_a", self.codomain_curve.a())
            .field("first_x_pullback", self.first.x_pullback())
            .field("first_y_pullback", self.first.y_pullback())
            .field("second_x_pullback", self.second.x_pullback())
            .field("second_y_pullback", self.second.y_pullback())
            .finish()
    }
}

/// Returns one valid short-Weierstrass pullback map.
///
/// The generated cases currently mix:
///
/// - self-maps on one curve given by identity or negation
/// - constant maps to rational finite points on an arbitrary codomain curve
pub fn arb_short_weierstrass_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<FunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    prop_oneof![
        arb_same_curve_function_field_map_case(curve_config),
        arb_constant_function_field_map_case(curve_config),
    ]
    .boxed()
}

/// Returns two composable pullback maps and their composition.
pub fn arb_composable_short_weierstrass_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<ComposableFunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    prop_oneof![
        arb_same_curve_composable_function_field_map_case(curve_config),
        arb_constant_chain_composable_function_field_map_case(curve_config),
    ]
    .boxed()
}

fn arb_same_curve_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<FunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(|curve| (Just(curve), any::<bool>()))
        .prop_map(|(curve, negate_y)| {
            let field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(curve.clone());
            let map = build_same_curve_map(curve.clone(), negate_y);

            FunctionFieldMapCase {
                domain_curve: curve.clone(),
                codomain_curve: curve.clone(),
                domain_field: crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(curve.clone()),
                codomain_field: field,
                map,
            }
        })
        .boxed()
}

fn arb_constant_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<FunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(move |domain_curve| {
            arb_curve_and_finite_point::<F>(curve_config).prop_map(
                move |(codomain_curve, point)| (domain_curve.clone(), codomain_curve, point),
            )
        })
        .prop_map(|(domain_curve, codomain_curve, point)| {
            let domain_field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(domain_curve.clone());
            let codomain_field = crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<F>::new(codomain_curve.clone());
            let map = build_constant_map(domain_curve.clone(), codomain_curve.clone(), &point);

            FunctionFieldMapCase {
                domain_curve,
                codomain_curve,
                domain_field,
                codomain_field,
                map,
            }
        })
        .boxed()
}

fn arb_same_curve_composable_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<ComposableFunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(|curve| (Just(curve), any::<bool>(), any::<bool>()))
        .prop_map(|(curve, first_negate, second_negate)| {
            let first = build_same_curve_map(curve.clone(), first_negate);
            let second = build_same_curve_map(curve.clone(), second_negate);
            let composite = first
                .compose(&second)
                .expect("same-curve self maps should compose");

            ComposableFunctionFieldMapCase {
                domain_curve: curve.clone(),
                middle_curve: curve.clone(),
                codomain_curve: curve,
                first,
                second,
                composite,
            }
        })
        .boxed()
}

fn arb_constant_chain_composable_function_field_map_case<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<ComposableFunctionFieldMapCase<F>>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(move |domain_curve| {
            arb_curve_and_finite_point::<F>(curve_config).prop_flat_map(
                move |(middle_curve, middle_point)| {
                    let domain_curve = domain_curve.clone();
                    arb_curve_and_finite_point::<F>(curve_config).prop_map(
                        move |(codomain_curve, codomain_point)| {
                            (
                                domain_curve.clone(),
                                middle_curve.clone(),
                                middle_point.clone(),
                                codomain_curve,
                                codomain_point,
                            )
                        },
                    )
                },
            )
        })
        .prop_map(
            |(domain_curve, middle_curve, middle_point, codomain_curve, codomain_point)| {
                let first =
                    build_constant_map(domain_curve.clone(), middle_curve.clone(), &middle_point);
                let second = build_constant_map(
                    middle_curve.clone(),
                    codomain_curve.clone(),
                    &codomain_point,
                );
                let composite = first
                    .compose(&second)
                    .expect("constant pullback maps should compose");

                ComposableFunctionFieldMapCase {
                    domain_curve,
                    middle_curve,
                    codomain_curve,
                    first,
                    second,
                    composite,
                }
            },
        )
        .boxed()
}

fn arb_curve_and_finite_point<F>(
    curve_config: CurveStrategyConfig,
) -> BoxedStrategy<(ShortWeierstrassCurve<F>, AffinePoint<F>)>
where
    F: EnumerableFiniteField + SqrtField + 'static,
    F::Elem: PartialEq + 'static,
    ShortWeierstrassCurve<F>:
        CurveModel<BaseField = F, Elem = F::Elem, Point = AffinePoint<F>> + EnumerableCurveModel,
{
    arb_nonsingular_curve::<F>(curve_config)
        .prop_flat_map(|curve| {
            let finite_points = curve
                .points()
                .into_iter()
                .filter(|point| !point.is_identity())
                .collect::<Vec<_>>();
            (Just(curve), prop::sample::select(finite_points))
        })
        .boxed()
}

fn build_same_curve_map<F>(
    curve: ShortWeierstrassCurve<F>,
    negate_y: bool,
) -> ShortWeierstrassFunctionFieldMap<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    let field =
        crate::elliptic_curves::short_weierstrass::function_fields::ShortWeierstrassFunctionField::<
            F,
        >::new(curve.clone());
    let y_pullback = if negate_y { field.y().neg() } else { field.y() };

    ShortWeierstrassFunctionFieldMap::new(curve.clone(), curve, field.x(), y_pullback)
        .expect("identity and negation should define valid self pullbacks")
}

fn build_constant_map<F>(
    domain_curve: ShortWeierstrassCurve<F>,
    codomain_curve: ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
) -> ShortWeierstrassFunctionFieldMap<F>
where
    F: Field,
    F::Elem: PartialEq,
{
    let AffinePoint::Finite { x, y } = point else {
        panic!("constant function-field maps use finite codomain points only");
    };

    let x_pullback = ShortWeierstrassFunction::<F>::from_rational_function(
        domain_curve.clone(),
        RationalFunction::<F>::constant(x.clone()),
    );
    let y_pullback = ShortWeierstrassFunction::<F>::from_rational_function(
        domain_curve.clone(),
        RationalFunction::<F>::constant(y.clone()),
    );

    ShortWeierstrassFunctionFieldMap::new(domain_curve, codomain_curve, x_pullback, y_pullback)
        .expect("finite codomain points should satisfy the codomain equation")
}

pub(crate) fn touch_function_field_map_case_fields() {
    let _ = |case: FunctionFieldMapCase<crate::fields::Fp17>| {
        let _ = case.domain_curve;
        let _ = case.codomain_curve;
        let _ = case.domain_field;
        let _ = case.codomain_field;
        let _ = case.map;
    };
    let _ = |case: ComposableFunctionFieldMapCase<crate::fields::Fp17>| {
        let _ = case.domain_curve;
        let _ = case.middle_curve;
        let _ = case.codomain_curve;
        let _ = case.first;
        let _ = case.second;
        let _ = case.composite;
    };
}
