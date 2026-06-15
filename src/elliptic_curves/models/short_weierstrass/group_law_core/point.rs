/// Internal affine point shape used by the shared short-Weierstrass formulas.
///
/// This deliberately mirrors the two public point surfaces that currently use
/// the same group law:
///
/// - `AffinePoint<F>` over the base field
/// - `ShortWeierstrassFunctionFieldPoint<F>` over `F(E)`
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ShortWeierstrassFormulaPoint<T> {
    Infinity,
    Affine { x: T, y: T },
}
