/// Strategy used to compute `E(ℚ)_tors`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RationalTorsionStrategy {
    /// Integral-model Lutz-Nagell candidate enumeration followed by exact
    /// Mazur-order verification.
    LutzNagell,
    /// Good reduction modulo a small prime, division-polynomial `x`-criteria,
    /// Hensel lifting of integer `x`-roots, exact `y`-recovery, and exact
    /// Mazur-order verification.
    GoodReductionHensel,
}
