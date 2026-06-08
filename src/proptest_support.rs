use core::fmt;

use proptest::prelude::*;

/// Shared property-testing strategies and small educational case builders.
///
/// Current scope:
/// - small prime-field elements and polynomial shapes
/// - small non-singular short-Weierstrass curves together with sampled
///   rational points
/// - one reusable quadratic extension and one small extension tower over `F17`
/// - tiny Vélu and composition cases over the fixed `F41` sample curve
///
/// Future work:
/// - widen the tower coverage beyond the current `F17 -> F17(sqrt(3)) -> degree-3`
///   path, ideally with multiple tower shapes and more shrink-friendly
///   semantic generators
/// - add composable Vélu cases over more than one small base curve instead of
///   concentrating the whole higher-level isogeny property surface on the
///   current `F41` family
/// - bias future isogeny strategies more deliberately toward special kernels
///   such as low-order, higher-codimension, or bridge-rich cases when those
///   become part of the educational surface
use crate::{
    elliptic_curves::{
        AffinePoint, CurveIsomorphism, CurveModel, EnumerableCurveModel, FiniteGroupCurveModel,
        GroupCurveModel, ShortWeierstrassCurve, ShortWeierstrassIsomorphism,
    },
    fields::{
        ExtensionField, ExtensionFieldElement, ExtensionFieldSpec, Field, Fp, FpElem,
        PolynomialModulus,
    },
    isogenies::{Isogeny, VeluIsogeny, maps_equal_exhaustively},
    polynomials::{
        DensePolynomial, Monomial, MultivariatePolynomial, MultivariateTerm, SparsePolynomial,
        SparsePolynomialTerm,
    },
};

type F17 = Fp<17>;
type F41 = Fp<41>;
type Curve41 = ShortWeierstrassCurve<F41>;
type Point41 = AffinePoint<F41>;

pub(crate) fn fp_elem<const P: u64>() -> impl Strategy<Value = FpElem<P>> {
    (0..P).prop_map(Fp::<P>::elem_from_u64)
}

pub(crate) fn nonzero_fp_elem<const P: u64>() -> impl Strategy<Value = FpElem<P>> {
    (1..P).prop_map(Fp::<P>::elem_from_u64)
}

pub(crate) fn dense_polynomial<const P: u64>(
    max_len: usize,
) -> impl Strategy<Value = DensePolynomial<Fp<P>>> {
    prop::collection::vec(0..P, 0..=max_len).prop_map(|coefficients| {
        DensePolynomial::<Fp<P>>::new(
            coefficients
                .into_iter()
                .map(Fp::<P>::elem_from_u64)
                .collect(),
        )
    })
}

pub(crate) fn sparse_polynomial<const P: u64>(
    max_terms: usize,
    max_degree: usize,
) -> impl Strategy<Value = SparsePolynomial<Fp<P>>> {
    prop::collection::vec((0..P, 0usize..=max_degree), 0..=max_terms).prop_map(|terms| {
        SparsePolynomial::<Fp<P>>::new(
            terms
                .into_iter()
                .map(|(coefficient, degree)| SparsePolynomialTerm {
                    coefficient: Fp::<P>::elem_from_u64(coefficient),
                    degree,
                })
                .collect(),
        )
    })
}

pub(crate) fn multivariate_polynomial<const P: u64>(
    arity: usize,
    max_terms: usize,
    max_exponent: usize,
) -> impl Strategy<Value = MultivariatePolynomial<Fp<P>>> {
    prop::collection::vec(
        (0..P, prop::collection::vec(0usize..=max_exponent, arity)),
        0..=max_terms,
    )
    .prop_map(move |terms| {
        MultivariatePolynomial::<Fp<P>>::new(
            arity,
            terms
                .into_iter()
                .map(|(coefficient, exponents)| MultivariateTerm {
                    coefficient: Fp::<P>::elem_from_u64(coefficient),
                    monomial: Monomial::new(exponents),
                })
                .collect(),
        )
        .expect("arity matches by construction")
    })
}

pub(crate) fn non_singular_short_weierstrass_curve<const P: u64>()
-> impl Strategy<Value = ShortWeierstrassCurve<Fp<P>>> {
    (0..P, 0..P).prop_filter_map("curve must be non-singular", |(a, b)| {
        ShortWeierstrassCurve::<Fp<P>>::new(Fp::<P>::elem_from_u64(a), Fp::<P>::elem_from_u64(b))
            .ok()
    })
}

pub(crate) fn curve_and_rational_point<const P: u64>()
-> BoxedStrategy<(ShortWeierstrassCurve<Fp<P>>, AffinePoint<Fp<P>>)> {
    non_singular_short_weierstrass_curve::<P>()
        .prop_flat_map(|curve| {
            let points = curve.points();
            let point_count = points.len();
            (Just(curve), Just(points), 0..point_count)
                .prop_map(|(curve, points, index)| (curve, points[index].clone()))
        })
        .boxed()
}

pub(crate) fn distinct_fp_elements<const P: u64>(count: usize) -> BoxedStrategy<Vec<FpElem<P>>> {
    prop::sample::subsequence(
        (0..P).map(Fp::<P>::elem_from_u64).collect::<Vec<_>>(),
        count..=count,
    )
    .boxed()
}

pub(crate) struct ProptestF17Sqrt3Spec;

impl ExtensionFieldSpec for ProptestF17Sqrt3Spec {
    type Base = F17;

    const NAME: &'static str = "proptest F17(sqrt(3))";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<Self::Base>::new(vec![
            <Self::Base as Field>::from_i64(-3),
            <Self::Base as Field>::zero(),
            <Self::Base as Field>::one(),
        ])
        .expect("x^2 - 3 should be a valid structural modulus")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Self::defining_modulus().check_field_modulus_requirements()
    }
}

pub(crate) type ProptestF17Sqrt3Field = ExtensionField<ProptestF17Sqrt3Spec>;

pub(crate) struct ProptestF17TowerSpec;

impl ExtensionFieldSpec for ProptestF17TowerSpec {
    type Base = ProptestF17Sqrt3Field;

    const NAME: &'static str = "proptest F17(sqrt(3))(u)";

    fn defining_modulus() -> PolynomialModulus<Self::Base> {
        PolynomialModulus::<ProptestF17Sqrt3Field>::new(vec![
            ProptestF17Sqrt3Field::one(),
            ProptestF17Sqrt3Field::one(),
            ProptestF17Sqrt3Field::zero(),
            ProptestF17Sqrt3Field::one(),
        ])
        .expect("tower modulus should be structurally valid")
    }

    fn check_field_conditions() -> Result<(), crate::fields::FieldError> {
        Ok(())
    }
}

pub(crate) type ProptestF17TowerField = ExtensionField<ProptestF17TowerSpec>;
pub(crate) type ProptestF17Sqrt3Elem = ExtensionFieldElement<ProptestF17Sqrt3Spec>;
pub(crate) type ProptestF17TowerElem = ExtensionFieldElement<ProptestF17TowerSpec>;

#[derive(Clone, Debug)]
pub(crate) struct TowerElementCase {
    pub base_left: FpElem<17>,
    pub base_right: FpElem<17>,
    pub quadratic_left: ProptestF17Sqrt3Elem,
    pub quadratic_right: ProptestF17Sqrt3Elem,
    pub tower_left: ProptestF17TowerElem,
    pub tower_right: ProptestF17TowerElem,
}

#[derive(Clone)]
pub(crate) struct CyclicKernelCase {
    pub curve: Curve41,
    pub generator: Point41,
    pub kernel_point: Point41,
    pub sample_point: Point41,
    pub coset_point: Point41,
    pub isogeny: VeluIsogeny<Curve41>,
}

impl fmt::Debug for CyclicKernelCase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CyclicKernelCase")
            .field("generator", &self.generator)
            .field("kernel_point", &self.kernel_point)
            .field("sample_point", &self.sample_point)
            .field("coset_point", &self.coset_point)
            .field("kernel_size", &self.isogeny.kernel_points().len())
            .finish()
    }
}

#[derive(Clone)]
pub(crate) struct ComposableVeluCase {
    pub first: VeluIsogeny<Curve41>,
    pub second_strict: VeluIsogeny<Curve41>,
    pub bridge: ShortWeierstrassIsomorphism<F41>,
    pub second_bridged: VeluIsogeny<Curve41>,
    pub sample_point: Point41,
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

pub(crate) fn tower_element_case() -> BoxedStrategy<TowerElementCase> {
    (
        fp_elem::<17>(),
        fp_elem::<17>(),
        semantic_quadratic_element(),
        semantic_quadratic_element(),
        semantic_tower_element(),
        semantic_tower_element(),
    )
        .prop_map(
            |(base_left, base_right, quadratic_left, quadratic_right, tower_left, tower_right)| {
                TowerElementCase {
                    base_left,
                    base_right,
                    quadratic_left,
                    quadratic_right,
                    tower_left,
                    tower_right,
                }
            },
        )
        .boxed()
}

pub(crate) fn cyclic_kernel_case() -> BoxedStrategy<CyclicKernelCase> {
    prop::sample::select(all_cyclic_kernel_cases()).boxed()
}

pub(crate) fn composable_velu_case() -> BoxedStrategy<ComposableVeluCase> {
    prop::sample::select(all_composable_velu_cases()).boxed()
}

fn semantic_quadratic_element() -> BoxedStrategy<ProptestF17Sqrt3Elem> {
    prop_oneof![
        Just(ProptestF17Sqrt3Field::zero()),
        fp_elem::<17>().prop_map(ProptestF17Sqrt3Field::from_base),
        nonzero_fp_elem::<17>()
            .prop_map(|coefficient| ProptestF17Sqrt3Field::element(vec![F17::zero(), coefficient])),
        (fp_elem::<17>(), nonzero_fp_elem::<17>())
            .prop_map(|(left, right)| ProptestF17Sqrt3Field::element(vec![left, right])),
    ]
    .boxed()
}

fn semantic_tower_element() -> BoxedStrategy<ProptestF17TowerElem> {
    prop_oneof![
        Just(ProptestF17TowerField::zero()),
        semantic_quadratic_element().prop_map(ProptestF17TowerField::from_base),
        semantic_quadratic_element()
            .prop_filter(
                "degree-one tower coefficient should be non-zero",
                |coefficient| { !ProptestF17Sqrt3Field::is_zero(coefficient) }
            )
            .prop_map(|coefficient| ProptestF17TowerField::element(vec![
                ProptestF17Sqrt3Field::zero(),
                coefficient
            ])),
        (
            semantic_quadratic_element(),
            semantic_quadratic_element().prop_filter(
                "tower element should have a non-zero higher coefficient",
                |coefficient| !ProptestF17Sqrt3Field::is_zero(coefficient),
            ),
            semantic_quadratic_element(),
        )
            .prop_map(|(constant, linear, quadratic)| {
                ProptestF17TowerField::element(vec![constant, linear, quadratic])
            }),
    ]
    .boxed()
}

fn proptest_curve41() -> Curve41 {
    ShortWeierstrassCurve::<F41>::new(F41::from_i64(2), F41::from_i64(3))
        .expect("shared F41 curve should stay valid")
}

fn all_cyclic_kernel_cases() -> Vec<CyclicKernelCase> {
    let curve = proptest_curve41();
    let all_points = curve.points();
    let mut kernels = Vec::<Vec<Point41>>::new();
    let mut cases = Vec::new();

    for generator in &all_points {
        if curve.is_identity(generator) {
            continue;
        }

        let Some(order) = curve.point_order(generator) else {
            continue;
        };
        if order <= 1 {
            continue;
        }

        let Ok(isogeny) = VeluIsogeny::from_generator(curve.clone(), generator.clone()) else {
            continue;
        };
        let kernel_points = isogeny.kernel_points().to_vec();
        if kernel_points.len() == all_points.len() {
            continue;
        }
        if kernels
            .iter()
            .any(|existing| same_point_set(existing, &kernel_points))
        {
            continue;
        }
        kernels.push(kernel_points.clone());

        let kernel_point = kernel_points
            .iter()
            .find(|point| !curve.is_identity(point))
            .expect("non-trivial cyclic kernel should contain a non-identity point")
            .clone();
        let sample_point = all_points
            .iter()
            .find(|point| !kernel_points.contains(point))
            .expect("sample curve should have points outside the kernel")
            .clone();
        let coset_point = curve
            .add(&sample_point, &kernel_point)
            .expect("kernel translation should stay on the curve");

        cases.push(CyclicKernelCase {
            curve: curve.clone(),
            generator: generator.clone(),
            kernel_point,
            sample_point,
            coset_point,
            isogeny,
        });
    }

    cases
}

fn all_composable_velu_cases() -> Vec<ComposableVeluCase> {
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

        for bridge in nontrivial_bridges(first_case.isogeny.codomain()) {
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

fn unique_velu_isogenies_on(curve: &Curve41) -> Vec<VeluIsogeny<Curve41>> {
    let mut kernels = Vec::<Vec<Point41>>::new();
    let mut isogenies = Vec::new();

    for point in curve.points() {
        if curve.is_identity(&point) {
            continue;
        }

        let Some(order) = curve.point_order(&point) else {
            continue;
        };
        if order <= 1 {
            continue;
        }

        let Ok(isogeny) = VeluIsogeny::from_generator(curve.clone(), point.clone()) else {
            continue;
        };
        let kernel = isogeny.kernel_points().to_vec();

        if kernels
            .iter()
            .any(|existing| same_point_set(existing, &kernel))
        {
            continue;
        }

        kernels.push(kernel);
        isogenies.push(isogeny);
    }

    isogenies
}

fn nontrivial_bridges(curve: &Curve41) -> Vec<ShortWeierstrassIsomorphism<F41>> {
    [F41::from_i64(2), F41::from_i64(3), F41::from_i64(5)]
        .into_iter()
        .filter_map(|u| ShortWeierstrassIsomorphism::new(curve.clone(), u).ok())
        .collect()
}

fn deduplicate_composable_cases(cases: Vec<ComposableVeluCase>) -> Vec<ComposableVeluCase> {
    let mut unique = Vec::new();

    for candidate in cases {
        let already_present = unique.iter().any(|existing: &ComposableVeluCase| {
            same_point_set(
                candidate.first.kernel_points(),
                existing.first.kernel_points(),
            ) && same_point_set(
                candidate.second_bridged.kernel_points(),
                existing.second_bridged.kernel_points(),
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

fn same_point_set(left: &[Point41], right: &[Point41]) -> bool {
    left.len() == right.len()
        && left.iter().all(|point| right.contains(point))
        && right.iter().all(|point| left.contains(point))
}
