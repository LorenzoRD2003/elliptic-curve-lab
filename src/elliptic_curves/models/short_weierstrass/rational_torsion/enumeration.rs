/// Bound used by the first rational-torsion route.
///
/// For `E(Q)_tors`, Mazur's theorem bounds every point order by the least
/// common multiple of the permitted cyclic orders:
/// `lcm(1, 2, ..., 10, 12) = 27720`.
pub(crate) const MAZUR_TORSION_EXPONENT_BOUND: usize = 27_720;

/// Cyclic orders permitted by Mazur's theorem over `Q`.
pub(crate) const MAZUR_CYCLIC_ORDERS: &[usize] = &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];

/// Product-family parameters `m` for `ℤ/2ℤ × ℤ/2mℤ`.
pub(crate) const MAZUR_PRODUCT_PARAMETERS: &[usize] = &[1, 2, 3, 4];
