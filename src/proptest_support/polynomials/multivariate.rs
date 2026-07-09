use crate::fields::traits::*;
use core::fmt::Debug;

use proptest::prelude::*;

use crate::fields::traits::EnumerableFiniteField;
use crate::polynomials::MultivariatePolynomial;
use crate::polynomials::multivariate::{Monomial, MultivariateTerm};
use crate::proptest_support::config::PolynomialStrategyConfig;
use crate::proptest_support::polynomials::shared::arb_field_elem;

/// Returns a multivariate polynomial over a small enumerable finite field.
pub fn arb_multivariate_polynomial<F>(
    config: PolynomialStrategyConfig,
) -> BoxedStrategy<MultivariatePolynomial<F>>
where
    F: EnumerableFiniteField + Debug + 'static,
{
    let arity = config.arity;
    prop::collection::vec(
        (
            arb_field_elem::<F>(),
            prop::collection::vec(0usize..=config.max_exponent, arity),
        ),
        0..=config.max_terms,
    )
    .prop_map(move |terms| {
        MultivariatePolynomial::<F>::new(
            arity,
            terms
                .into_iter()
                .map(|(coefficient, exponents)| {
                    MultivariateTerm::new(coefficient, Monomial::new(exponents))
                })
                .collect(),
        )
        .expect("arity matches by construction")
    })
    .boxed()
}
