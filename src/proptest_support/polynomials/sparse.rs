use crate::fields::traits::*;
use core::fmt::Debug;

use proptest::prelude::*;

use crate::fields::traits::EnumerableFiniteField;
use crate::polynomials::SparsePolynomial;
use crate::polynomials::sparse::SparsePolynomialTerm;
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::shared::arb_field_elem;

/// Returns a sparse polynomial over a small enumerable finite field.
pub fn arb_sparse_polynomial<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<SparsePolynomial<F>>
where
    F: EnumerableFiniteField + Debug + 'static,
{
    prop::collection::vec(
        (arb_field_elem::<F>(), 0usize..=config.max_degree),
        0..=config.max_terms,
    )
    .prop_map(|terms| {
        SparsePolynomial::<F>::new(
            terms
                .into_iter()
                .map(|(coefficient, degree)| SparsePolynomialTerm {
                    coefficient,
                    degree,
                })
                .collect(),
        )
    })
    .boxed()
}
