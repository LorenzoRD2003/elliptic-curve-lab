use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::short_weierstrass::ShortWeierstrassCurve;
use crate::elliptic_curves::traits::CurveModel;
use crate::fields::FiniteField;

/// Returns the `p^k`-power Frobenius twist of a short-Weierstrass curve.
///
/// If `E : y^2 = x^3 + ax + b` over a finite field of characteristic `p`,
/// this helper returns the `p^k`-power Frobenius twist
/// `E^(p^k) : y^2 = x^3 + a^(p^k)x + b^(p^k)`.
///
/// When `a` and `b` already lie in `F_p`, this twist agrees with the original
/// curve. Over a larger field such as `F_{p^r}`, the coefficients may move.
///
/// The implementation reduces `k` modulo the extension degree of the chosen
/// finite field backend. If `F` represents `F_{p^n}`, then `π_p^n` is already
/// the identity on every field element, so the coefficient action depends only
/// on `k mod n`.
pub fn frobenius_twist_power<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    power: u32,
) -> Result<ShortWeierstrassCurve<F>, CurveError> {
    ShortWeierstrassCurve::<F>::new(
        frobenius_power_element::<F>(curve.a(), power),
        frobenius_power_element::<F>(curve.b(), power),
    )
}

/// Applies the absolute Frobenius `π_p^k` to a point's coordinates.
///
/// The input point must lie on `curve`. The returned coordinates lie on the
/// `p^k`-power Frobenius twist returned by [`frobenius_twist_power`].
///
/// Concretely, this helper applies the coordinate transformation
///
/// - `π_p^k(O) = O`
/// - `π_p^k(x, y) = (x^(p^k), y^(p^k))`
///
/// So if `P = (x, y)` lies on `E : y^2 = x^3 + ax + b`,
/// then the returned point lies on
/// `E^(p^k) : y^2 = x^3 + a^(p^k)x + b^(p^k)`.
///
/// As with [`frobenius_twist_power`], the implementation reduces `k` modulo
/// the extension degree of the represented finite field, since over
/// `F_{p^n}` one already has `π_p^n = id` on field elements.
pub fn absolute_frobenius_power_point<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    power: u32,
) -> Result<AffinePoint<F>, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }
    Ok(point.map_coordinates::<F, _>(|coordinate| frobenius_power_element::<F>(coordinate, power)))
}

/// Applies the relative Frobenius `π_q` to a point on a curve over `F_q`.
///
/// When `F` has size `q = p^r`, the relative Frobenius acts by the coordinate
/// transformation
///
/// - `π_q(O) = O`
/// - `π_q(x, y) = (x^q, y^q)`
///
/// and returns another point on the same short-Weierstrass curve.
///
/// This is the curve-side map whose fixed points are meant to model
/// `E(F_q)` inside the currently represented finite base field.
///
/// Since the curve is assumed to be defined over `F_q`, this relative
/// Frobenius lands back on the same curve rather than on a distinct
/// Frobenius twist.
///
/// For the current backend model, where `F` itself is the represented field
/// `F_q`, this means the coordinate action is already the identity on every
/// field element. The implementation therefore returns the input point after
/// membership validation.
///
/// In particular, this public helper does not take an iterate parameter:
/// once the ambient backend already is `F_q`, the observable point action of
/// `π_q` is the same identity map that every power `π_q^k` would induce.
///
/// Complexity: Θ(1)
pub fn relative_frobenius_point<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
) -> Result<AffinePoint<F>, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    Ok(point.clone())
}

fn frobenius_power_element<F: FiniteField>(element: &F::Elem, power: u32) -> F::Elem {
    let extension_degree = F::extension_degree().get();
    let reduced_power = power % extension_degree;
    let mut image = element.clone();
    for _ in 0..reduced_power {
        image = F::pow(&image, F::characteristic());
    }
    image
}
