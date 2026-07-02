use num_bigint::BigUint;

use crate::elliptic_curves::{
    group_algorithms::cyclic_roots::{
        CyclicPrimeRootError, CyclicPrimeRootReport, algorithm::compute_cyclic_prime_root_report,
    },
    traits::GroupCurveModel,
};

/// Curve models that support the staged cyclic prime-root algorithm.
///
/// The group law is written additively. For a point `γ`, this solves
/// `[r]ρ = γ` in a finite cyclic curve group of known order `|G|`, using a
/// caller-supplied generator `δ` for the `r`-Sylow subgroup.
pub(crate) trait CyclicPrimeRootCurveModel: GroupCurveModel
where
    Self::Point: Clone + PartialEq,
{
    /// Attempts to compute one `r`-th root of `target`.
    fn cyclic_prime_root(
        &self,
        target: &Self::Point,
        root_degree: BigUint,
        group_order: BigUint,
        sylow_generator: &Self::Point,
    ) -> Result<CyclicPrimeRootReport<Self::Point>, CyclicPrimeRootError> {
        compute_cyclic_prime_root_report(self, target, root_degree, group_order, sylow_generator)
    }
}

impl<C: GroupCurveModel> CyclicPrimeRootCurveModel for C where C::Point: Clone + PartialEq {}
