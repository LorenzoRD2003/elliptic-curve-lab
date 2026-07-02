use num_bigint::BigUint;

use crate::elliptic_curves::{
    group_algorithms::cyclic_roots::{
        CyclicPrimeRootError, CyclicPrimeRootReport, algorithm::compute_cyclic_prime_root_report,
    },
    traits::GroupCurveModel,
};

/// Curve models used as externally certified finite cyclic groups.
///
/// The group law is written additively. For a point `γ`, this solves
/// `[r]ρ = γ` in a finite cyclic curve group of known order `|G|`, using a
/// caller-supplied generator `δ` for the `r`-Sylow subgroup.
///
/// The exercise relies on the cyclic-group order identity
/// `|nγ| = |γ| / gcd(n, |γ|)`: after writing `|G| = a r^k`, the point
/// `α = aγ` lands in the `r`-Sylow subgroup, so a generator `δ` of that
/// subgroup lets the algorithm search for `x` with `α = xδ`.
///
/// Important: this trait does not prove that the represented curve group is
/// cyclic. It is a capability surface for callers that already know, from
/// external data or a separate certificate, that the finite curve group is
/// cyclic and that the supplied `δ` generates the full `r`-Sylow subgroup.
/// On a non-cyclic group the algorithm can return `NoRoot` or a missing-log
/// error even when roots exist outside `<δ>`.
pub trait CyclicGroupPrimeRootCurveModel: GroupCurveModel
where
    Self::Point: Clone + PartialEq,
{
    /// Attempts to compute one `r`-th root of `target`.
    ///
    /// The caller is responsible for the cyclic-group hypothesis. In
    /// particular, when `r | |G|`, `sylow_generator` must generate the full
    /// `r`-Sylow subgroup, not merely have some `r`-power order.
    ///
    /// In the nontrivial Sylow case, a successful route returns
    /// `ρ = s(x/r)δ + tβ`, where `β = r^kγ` and
    /// `s a + t r^(k+1) = 1`.
    fn cyclic_group_prime_root(
        &self,
        target: &Self::Point,
        root_degree: BigUint,
        group_order: BigUint,
        sylow_generator: &Self::Point,
    ) -> Result<CyclicPrimeRootReport<Self::Point>, CyclicPrimeRootError> {
        compute_cyclic_prime_root_report(self, target, root_degree, group_order, sylow_generator)
    }
}

impl<C: GroupCurveModel> CyclicGroupPrimeRootCurveModel for C where C::Point: Clone + PartialEq {}
