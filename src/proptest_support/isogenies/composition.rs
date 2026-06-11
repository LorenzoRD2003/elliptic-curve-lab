use core::fmt;

use proptest::prelude::*;

use crate::elliptic_curves::{
    CurveIsomorphism, EnumerableCurveModel, FiniteGroupCurveModel, ShortWeierstrassIsomorphism,
};
use crate::fields::Field;
use crate::isogenies::{Isogeny, VeluIsogeny, maps_equal_exhaustively};
use crate::proptest_support::combinators::same_membership_set;
use crate::proptest_support::config::IsogenyStrategyConfig;
use crate::proptest_support::isogenies::kernels::{Curve41, F41, unique_velu_isogenies_on};
use crate::proptest_support::isogenies::velu::all_cyclic_kernel_cases;

/// Composable tiny Vélu case, optionally bridged by a short-Weierstrass
/// scaling isomorphism.
#[derive(Clone)]
pub struct ComposableVeluCase {
    pub first: VeluIsogeny<Curve41>,
    pub second_strict: VeluIsogeny<Curve41>,
    pub bridge: ShortWeierstrassIsomorphism<F41>,
    pub second_bridged: VeluIsogeny<Curve41>,
    pub sample_point: crate::elliptic_curves::AffinePoint<F41>,
}

impl fmt::Debug for ComposableVeluCase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ComposableVeluCase")
            .field("first_kernel_size", &self.first.kernel_points().len())
            .field(
                "second_strict_kernel_size",
                &self.second_strict.kernel_points().len(),
            )
            .field(
                "second_bridged_kernel_size",
                &self.second_bridged.kernel_points().len(),
            )
            .field("bridge_scale", self.bridge.scaling_factor())
            .field("sample_point", &self.sample_point)
            .finish()
    }
}

/// Returns one tiny composable Vélu fixture over the canonical `F41` family.
pub fn arb_composable_velu_case() -> BoxedStrategy<ComposableVeluCase> {
    prop::sample::select(all_composable_velu_cases(IsogenyStrategyConfig::default())).boxed()
}

pub(crate) fn all_composable_velu_cases(config: IsogenyStrategyConfig) -> Vec<ComposableVeluCase> {
    let mut cases = Vec::new();

    for first_case in all_cyclic_kernel_cases() {
        let sample_point = first_case.sample_point.clone();
        let strict_second = unique_velu_isogenies_on(first_case.isogeny.codomain());

        for second_strict in strict_second {
            cases.push(ComposableVeluCase {
                first: first_case.isogeny.clone(),
                second_strict,
                bridge: ShortWeierstrassIsomorphism::new(
                    first_case.isogeny.codomain().clone(),
                    F41::one(),
                )
                .expect("identity scaling should build"),
                second_bridged: VeluIsogeny::from_generator(
                    first_case.isogeny.codomain().clone(),
                    first_case
                        .isogeny
                        .codomain()
                        .points_of_order(2)
                        .into_iter()
                        .next()
                        .or_else(|| {
                            first_case
                                .isogeny
                                .codomain()
                                .points()
                                .into_iter()
                                .find(|point| {
                                    first_case
                                        .isogeny
                                        .codomain()
                                        .point_order(point)
                                        .is_some_and(|order| order > 1)
                                })
                        })
                        .expect("codomain should contain a non-trivial torsion point"),
                )
                .expect("fallback bridged isogeny should build"),
                sample_point: sample_point.clone(),
            });
        }

        for bridge in nontrivial_bridges(first_case.isogeny.codomain(), config) {
            for second_bridged in unique_velu_isogenies_on(bridge.codomain()) {
                cases.push(ComposableVeluCase {
                    first: first_case.isogeny.clone(),
                    second_strict: VeluIsogeny::from_generator(
                        first_case.isogeny.codomain().clone(),
                        first_case
                            .isogeny
                            .codomain()
                            .points()
                            .into_iter()
                            .find(|point| {
                                first_case
                                    .isogeny
                                    .codomain()
                                    .point_order(point)
                                    .is_some_and(|order| order > 1)
                            })
                            .expect("codomain should contain a non-trivial point"),
                    )
                    .expect("strict fallback isogeny should build"),
                    bridge: bridge.clone(),
                    second_bridged,
                    sample_point: sample_point.clone(),
                });
            }
        }
    }

    deduplicate_composable_cases(cases)
}

fn nontrivial_bridges(
    curve: &Curve41,
    config: IsogenyStrategyConfig,
) -> Vec<ShortWeierstrassIsomorphism<F41>> {
    config
        .preferred_bridge_scales
        .into_iter()
        .map(F41::from_i64)
        .filter_map(|u| ShortWeierstrassIsomorphism::new(curve.clone(), u).ok())
        .collect()
}

fn deduplicate_composable_cases(cases: Vec<ComposableVeluCase>) -> Vec<ComposableVeluCase> {
    let mut unique = Vec::new();

    for candidate in cases {
        let already_present = unique.iter().any(|existing: &ComposableVeluCase| {
            same_membership_set(
                &candidate.first.kernel_points(),
                &existing.first.kernel_points(),
            ) && same_membership_set(
                &candidate.second_bridged.kernel_points(),
                &existing.second_bridged.kernel_points(),
            ) && candidate.bridge.scaling_factor() == existing.bridge.scaling_factor()
                && maps_equal_exhaustively::<_, _, Curve41, Curve41>(
                    &candidate.second_strict,
                    &existing.second_strict,
                )
                .unwrap_or(false)
        });

        if !already_present {
            unique.push(candidate);
        }
    }

    unique
}
