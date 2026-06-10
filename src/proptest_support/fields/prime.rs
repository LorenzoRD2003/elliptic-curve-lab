use proptest::prelude::*;

use crate::fields::{Field, Fp, FpElem};

/// Returns a strategy for arbitrary elements of `GF(P)`.
pub fn arb_fp_elem<const P: u64>() -> BoxedStrategy<FpElem<P>> {
    (0..P).prop_map(Fp::<P>::elem_from_u64).boxed()
}

/// Returns a strategy for non-zero elements of `GF(P)`.
pub fn arb_nonzero_fp_elem<const P: u64>() -> BoxedStrategy<FpElem<P>> {
    (1..P).prop_map(Fp::<P>::elem_from_u64).boxed()
}

/// Returns a strategy for `count` pairwise distinct elements of `GF(P)`.
pub fn arb_distinct_fp_elems<const P: u64>(count: usize) -> BoxedStrategy<Vec<FpElem<P>>> {
    prop::sample::subsequence(
        (0..P).map(Fp::<P>::elem_from_u64).collect::<Vec<_>>(),
        count..=count,
    )
    .boxed()
}
