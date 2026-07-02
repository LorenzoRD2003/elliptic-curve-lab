use crate::fields::traits::*;
use core::fmt::Debug;

use proptest::prelude::*;

/// Returns a strategy that selects one element from a non-empty vector.
pub(crate) fn choose_from<T: Clone + Debug + 'static>(values: Vec<T>) -> BoxedStrategy<T> {
    prop::sample::select(values).boxed()
}

/// Returns `true` when both slices represent the same finite set by
/// membership, ignoring order.
pub(crate) fn same_membership_set<T: PartialEq>(left: &[T], right: &[T]) -> bool {
    left.len() == right.len()
        && left.iter().all(|value| right.contains(value))
        && right.iter().all(|value| left.contains(value))
}

pub(crate) fn touch_combinator_inventory() {
    let _ = choose_from(vec![1u8, 2u8, 3u8]);
    let _ = same_membership_set(&[1u8, 2u8], &[2u8, 1u8]);
}
