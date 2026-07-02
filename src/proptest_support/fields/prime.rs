use crate::fields::traits::*;
use proptest::prelude::*;

use crate::fields::traits::EnumerableFiniteField;

/// Returns a strategy for arbitrary elements of an enumerable prime field.
pub fn arb_fp_elem<F>() -> BoxedStrategy<F::Elem>
where
    F: EnumerableFiniteField,
    F::Elem: 'static,
{
    prop::sample::select(F::elements()).boxed()
}

/// Returns a strategy for non-zero elements of an enumerable prime field.
pub fn arb_nonzero_fp_elem<F>() -> BoxedStrategy<F::Elem>
where
    F: EnumerableFiniteField,
    F::Elem: 'static,
{
    prop::sample::select(
        F::elements()
            .into_iter()
            .filter(|element| !F::is_zero(element))
            .collect::<Vec<_>>(),
    )
    .boxed()
}

/// Returns a strategy for `count` pairwise distinct field elements.
pub fn arb_distinct_fp_elems<F>(count: usize) -> BoxedStrategy<Vec<F::Elem>>
where
    F: EnumerableFiniteField,
    F::Elem: 'static,
{
    prop::sample::subsequence(F::elements(), count..=count).boxed()
}
