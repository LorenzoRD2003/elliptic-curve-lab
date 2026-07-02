use crate::fields::traits::*;
use proptest::prelude::*;

use crate::fields::traits::EnumerableFiniteField;
use crate::proptest_support::combinators::choose_from;

/// Returns a field-element strategy by exhaustive enumeration over a small
/// enumerable finite field.
pub fn arb_field_elem<F>() -> BoxedStrategy<F::Elem>
where
    F: EnumerableFiniteField + 'static,
{
    choose_from(F::elements())
}

/// Returns a non-zero field-element strategy by exhaustive enumeration over a
/// small enumerable finite field.
pub fn arb_nonzero_field_elem<F>() -> BoxedStrategy<F::Elem>
where
    F: EnumerableFiniteField + 'static,
{
    choose_from(
        F::elements()
            .into_iter()
            .filter(|element| !F::is_zero(element))
            .collect(),
    )
}
