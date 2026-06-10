use proptest::prelude::*;

use crate::fields::EnumerableFiniteField;
use crate::proptest_support::combinators::choose_from;

/// Returns a field-element strategy by exhaustive enumeration over a small
/// enumerable finite field.
pub fn arb_field_elem<F>() -> BoxedStrategy<F::Elem>
where
    F: EnumerableFiniteField + 'static,
{
    choose_from(F::elements())
}
