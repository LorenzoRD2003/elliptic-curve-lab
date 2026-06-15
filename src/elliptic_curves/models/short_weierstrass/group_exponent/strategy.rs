use crate::elliptic_curves::short_weierstrass::point_order::PointOrderStrategy;

/// Public strategy choices for recovering or estimating `λ(E(F_q))`.
///
/// For a finite abelian group `G`, the exponent `λ(G) = lcm({|g| : g ∈ G})`
/// is also the maximum element order.
///
/// The current implementation distinguishes:
/// - [`Self::Exhaustive`], which computes the exact exponent on a tiny
///   enumerable curve group
/// - [`Self::RandomPoints`], which samples points with replacement, computes
///   their exact orders through one requested point-order strategy, and
///   accumulates the running least common multiple as a candidate for `λ(G)`
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupExponentStrategy {
    Exhaustive,
    RandomPoints {
        max_samples: usize,
        point_order_strategy: PointOrderStrategy,
    },
}
