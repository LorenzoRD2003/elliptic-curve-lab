use proptest::prelude::*;

use crate::fields::{
    EnumerableFiniteField, ExtensionField, ExtensionFieldElement, ExtensionFieldSpec, Field,
};
use crate::proptest_support::combinators::choose_from;

fn arb_base_elem<S>() -> BoxedStrategy<<S::Base as Field>::Elem>
where
    S: ExtensionFieldSpec + 'static,
    S::Base: EnumerableFiniteField,
{
    choose_from(S::Base::elements())
}

fn arb_nonzero_base_elem<S>() -> BoxedStrategy<<S::Base as Field>::Elem>
where
    S: ExtensionFieldSpec + 'static,
    S::Base: EnumerableFiniteField,
{
    arb_base_elem::<S>()
        .prop_filter("base coefficient should be non-zero", |value| {
            !S::Base::is_zero(value)
        })
        .boxed()
}

/// Returns a generic reduced representative in `Base[x] / (m(x))`.
pub fn arb_extension_elem<S>() -> BoxedStrategy<ExtensionFieldElement<S>>
where
    S: ExtensionFieldSpec + 'static,
    S::Base: EnumerableFiniteField,
{
    let degree = ExtensionField::<S>::extension_degree().get() as usize;
    prop::collection::vec(arb_base_elem::<S>(), 0..=degree)
        .prop_map(ExtensionField::<S>::element)
        .boxed()
}

/// Returns a shrink-friendly extension-field sample that prefers zero,
/// base embeddings, pure monomials, and only then denser representatives.
pub fn arb_semantic_extension_elem<S>() -> BoxedStrategy<ExtensionFieldElement<S>>
where
    S: ExtensionFieldSpec + 'static,
    S::Base: EnumerableFiniteField,
{
    let degree = ExtensionField::<S>::extension_degree().get() as usize;

    let zero = Just(ExtensionField::<S>::zero());
    let embedded_base = arb_base_elem::<S>().prop_map(ExtensionField::<S>::from_base);
    let dense = prop::collection::vec(arb_base_elem::<S>(), degree)
        .prop_filter(
            "dense representative should have a non-zero higher coefficient",
            |coefficients| {
                coefficients
                    .iter()
                    .enumerate()
                    .any(|(index, coefficient)| index > 0 && !S::Base::is_zero(coefficient))
            },
        )
        .prop_map(ExtensionField::<S>::element);

    if degree <= 1 {
        prop_oneof![zero, embedded_base].boxed()
    } else {
        let monomial =
            (1usize..degree, arb_nonzero_base_elem::<S>()).prop_map(move |(index, coefficient)| {
                let mut coefficients = vec![S::Base::zero(); degree];
                coefficients[index] = coefficient;
                ExtensionField::<S>::element(coefficients)
            });

        prop_oneof![zero, embedded_base, monomial, dense].boxed()
    }
}
