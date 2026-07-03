/// Bound used by the first rational-torsion route.
///
/// For `E(Q)_tors`, Mazur's theorem bounds every non-identity point order by
/// the least common multiple of the permitted nontrivial cyclic orders:
/// `lcm(2, ..., 10, 12) = 27720`.
pub(crate) const MAZUR_TORSION_EXPONENT_BOUND: usize = 27_720;

/// Nontrivial cyclic group orders permitted by Mazur's theorem over `Q`.
///
/// The identity has exact order `1`, but the order-one group is represented by
/// `RationalTorsionGroupShape::Trivial`, not by `Cyclic { order: 1 }`.
pub(crate) const MAZUR_CYCLIC_ORDERS: &[usize] = &[2, 3, 4, 5, 6, 7, 8, 9, 10, 12];

/// Product-family parameters `m` for `ℤ/2ℤ × ℤ/2mℤ`.
pub(crate) const MAZUR_PRODUCT_PARAMETERS: &[usize] = &[1, 2, 3, 4];
