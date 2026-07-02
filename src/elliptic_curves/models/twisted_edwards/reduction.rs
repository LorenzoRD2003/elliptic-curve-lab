use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::{AffinePoint, MontgomeryCurve, TwistedEdwardsCurve};

/// Errors for the classical birational transport between the affine
/// twisted-Edwards chart
///
/// `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2`
///
/// and the affine Montgomery chart
///
/// `E_{A,B}: B v^2 = u^3 + A u^2 + u`.
///
/// These are intentionally *point-level birational* errors, not whole-curve
/// conversion errors:
///
/// - `PointNotOnTwistedEdwards` / `PointNotOnMontgomery` mean the supplied
///   affine point is not on the claimed source curve at all
/// - `ExceptionalTwistedEdwardsPoint` / `ExceptionalMontgomeryPoint` mean the
///   point is on the source curve, but lies outside the affine open where the
///   classical rational formulas are defined
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TwistedEdwardsBirationalMapError {
    PointNotOnTwistedEdwards,
    PointNotOnMontgomery,
    ExceptionalTwistedEdwardsPoint,
    ExceptionalMontgomeryPoint,
}

impl fmt::Display for TwistedEdwardsBirationalMapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PointNotOnTwistedEdwards => {
                write!(
                    formatter,
                    "the supplied point does not lie on the twisted-Edwards source curve"
                )
            }
            Self::PointNotOnMontgomery => {
                write!(
                    formatter,
                    "the supplied point does not lie on the Montgomery source curve"
                )
            }
            Self::ExceptionalTwistedEdwardsPoint => {
                write!(
                    formatter,
                    "the supplied twisted-Edwards point lies outside the affine open where the birational map is defined"
                )
            }
            Self::ExceptionalMontgomeryPoint => {
                write!(
                    formatter,
                    "the supplied Montgomery point lies outside the affine open where the birational map is defined"
                )
            }
        }
    }
}

impl std::error::Error for TwistedEdwardsBirationalMapError {}

pub(crate) struct TwistedEdwardsMontgomeryReduction<F: Field> {
    twisted_edwards_curve: TwistedEdwardsCurve<F>,
    montgomery_curve: MontgomeryCurve<F>,
}

impl<F: Field> TwistedEdwardsMontgomeryReduction<F>
where
    F::Elem: Clone,
{
    pub(crate) fn from_twisted_edwards(twisted_edwards_curve: TwistedEdwardsCurve<F>) -> Self {
        let a_minus_d = F::sub(twisted_edwards_curve.a(), twisted_edwards_curve.d());
        let montgomery_curve = MontgomeryCurve::new(
            F::div(
                &F::mul(
                    &F::from_i64(2),
                    &F::add(twisted_edwards_curve.a(), twisted_edwards_curve.d()),
                ),
                &a_minus_d,
            )
            .expect("validated twisted-Edwards curve has a - d != 0"),
            F::div(&F::from_i64(4), &a_minus_d)
                .expect("validated twisted-Edwards curve has a - d != 0"),
        )
        .expect("the Montgomery companion of a validated twisted-Edwards curve is non-singular");

        Self {
            twisted_edwards_curve,
            montgomery_curve,
        }
    }

    pub(crate) fn from_montgomery(montgomery_curve: MontgomeryCurve<F>) -> Self {
        let twisted_edwards_curve = TwistedEdwardsCurve::new(
            F::div(
                &F::add(montgomery_curve.a(), &F::from_i64(2)),
                montgomery_curve.b(),
            )
            .expect("validated Montgomery curve has B != 0"),
            F::div(
                &F::sub(montgomery_curve.a(), &F::from_i64(2)),
                montgomery_curve.b(),
            )
            .expect("validated Montgomery curve has B != 0"),
        )
        .expect("the twisted-Edwards companion of a validated Montgomery curve is non-singular");

        Self {
            twisted_edwards_curve,
            montgomery_curve,
        }
    }

    fn transport_point_to_montgomery_open(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        if !self.twisted_edwards_curve.contains_affine_point(point) {
            return Err(TwistedEdwardsBirationalMapError::PointNotOnTwistedEdwards);
        }

        let AffinePoint::Finite { x, y } = point else {
            return Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint);
        };

        if F::is_zero(x) || F::eq(y, &F::one()) {
            return Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint);
        }

        let one_plus_y = F::add(&F::one(), y);
        let one_minus_y = F::sub(&F::one(), y);
        let montgomery_x = F::div(&one_plus_y, &one_minus_y)
            .expect("y != 1 on the twisted-Edwards birational open");
        let montgomery_y = F::div(&one_plus_y, &F::mul(x, &one_minus_y))
            .expect("x != 0 and y != 1 on the twisted-Edwards birational open");
        let image = AffinePoint::new(montgomery_x, montgomery_y);

        if self.montgomery_curve.contains_affine_point(&image) {
            Ok(image)
        } else {
            Err(TwistedEdwardsBirationalMapError::PointNotOnMontgomery)
        }
    }

    fn transport_point_to_montgomery_total(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        if !self.twisted_edwards_curve.contains_affine_point(point) {
            return Err(TwistedEdwardsBirationalMapError::PointNotOnTwistedEdwards);
        }

        let AffinePoint::Finite { x, y } = point else {
            return Err(TwistedEdwardsBirationalMapError::ExceptionalTwistedEdwardsPoint);
        };

        if F::is_zero(x) {
            let image = if F::eq(y, &F::one()) {
                AffinePoint::Infinity
            } else {
                AffinePoint::new(F::zero(), F::zero())
            };

            if self.montgomery_curve.contains_affine_point(&image) {
                return Ok(image);
            }

            return Err(TwistedEdwardsBirationalMapError::PointNotOnMontgomery);
        }

        self.transport_point_to_montgomery_open(point)
    }

    fn transport_point_from_montgomery_open(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        if !self.montgomery_curve.contains_affine_point(point) {
            return Err(TwistedEdwardsBirationalMapError::PointNotOnMontgomery);
        }

        let AffinePoint::Finite { x, y } = point else {
            return Err(TwistedEdwardsBirationalMapError::ExceptionalMontgomeryPoint);
        };

        if F::is_zero(y) || F::eq(x, &F::from_i64(-1)) {
            return Err(TwistedEdwardsBirationalMapError::ExceptionalMontgomeryPoint);
        }

        let twisted_x = F::div(x, y).expect("y != 0 on the Montgomery birational open");
        let twisted_y = F::div(&F::sub(x, &F::one()), &F::add(x, &F::one()))
            .expect("x != -1 on the Montgomery birational open");
        let image = AffinePoint::new(twisted_x, twisted_y);

        if self.twisted_edwards_curve.contains_affine_point(&image) {
            Ok(image)
        } else {
            Err(TwistedEdwardsBirationalMapError::PointNotOnTwistedEdwards)
        }
    }
}

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    /// Returns the canonically associated Montgomery companion
    ///
    /// `E_{A,B}: B y^2 = x^3 + A x^2 + x`
    ///
    /// with
    ///
    /// `A = 2(a + d) / (a - d)`, `B = 4 / (a - d)`.
    ///
    /// This is a *whole-curve* conversion, not a point conversion. The
    /// resulting Montgomery model is the standard companion of the same
    /// elliptic curve, and the formulas are defined globally because a
    /// validated twisted-Edwards descriptor satisfies `a ≠ d`.
    pub fn as_montgomery(&self) -> MontgomeryCurve<F> {
        MontgomeryCurve::from(self)
    }

    /// Transports one affine point through the classical birational map to the
    /// Montgomery companion on the affine open where the formulas are defined.
    ///
    /// This is intentionally *not* a total point-conversion API. It uses the
    /// rational formulas
    ///
    /// `u = (1 + y) / (1 - y)`, `v = (1 + y) / (x (1 - y))`,
    ///
    /// from the affine twisted-Edwards chart `(x, y)` to the affine
    /// Montgomery chart `(u, v)`. It is defined for `x ≠ 0`, `y ≠ 1`,
    ///
    /// so this method rejects exceptional twisted-Edwards points such as the
    /// neutral element `(0, 1)` and the second `x = 0` point `(0, -1)`.
    pub fn try_point_to_montgomery_open(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        TwistedEdwardsMontgomeryReduction::from_twisted_edwards(self.clone())
            .transport_point_to_montgomery_open(point)
    }

    /// Transports one affine twisted-Edwards point to the Montgomery
    /// companion, extending the classical birational formulas at the two
    /// exceptional `x = 0` points.
    ///
    /// This agrees with the birational-open map when `x ≠ 0` and `y ≠ 1`,
    /// and additionally sets
    ///
    /// - `(0, 1) -> O`, `(0, -1) -> (0, 0)`
    ///
    /// so this direction is total on affine twisted-Edwards points.
    pub fn point_to_montgomery(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        TwistedEdwardsMontgomeryReduction::from_twisted_edwards(self.clone())
            .transport_point_to_montgomery_total(point)
    }
}

impl<F: Field> MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    /// Returns the canonically associated twisted-Edwards companion
    ///
    /// `E_{a,d}: a x^2 + y^2 = 1 + d x^2 y^2`
    ///
    /// with `a = (A + 2) / B`, `d = (A - 2) / B`.
    ///
    /// This is a *whole-curve* conversion, not a point conversion. The
    /// formulas are defined globally because a validated Montgomery descriptor
    /// satisfies `B ≠ 0`.
    pub fn as_twisted_edwards(&self) -> TwistedEdwardsCurve<F> {
        TwistedEdwardsCurve::from(self)
    }

    /// Transports one affine point through the classical birational map to the
    /// twisted-Edwards companion on the affine open where the formulas are
    /// defined.
    ///
    /// This is intentionally *not* a total point-conversion API. It uses the
    /// inverse rational formulas
    ///
    /// `x = u / v`, `y = (u - 1) / (u + 1)`,
    ///
    /// from the affine Montgomery chart `(u, v)` to the affine twisted-Edwards
    /// chart `(x, y)`. It is defined for `v ≠ 0`, `u ≠ -1` and, because this API
    /// works only on the affine Montgomery chart, the point at infinity is also excluded.
    ///
    /// So this method rejects exceptional Montgomery points such as `O`,
    /// affine points with `y = 0`, and affine points with `x = -1`.
    pub fn try_point_to_twisted_edwards_open(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, TwistedEdwardsBirationalMapError> {
        TwistedEdwardsMontgomeryReduction::from_montgomery(self.clone())
            .transport_point_from_montgomery_open(point)
    }
}

impl<F: Field> From<&TwistedEdwardsCurve<F>> for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: &TwistedEdwardsCurve<F>) -> Self {
        TwistedEdwardsMontgomeryReduction::from_twisted_edwards(curve.clone()).montgomery_curve
    }
}

impl<F: Field> From<TwistedEdwardsCurve<F>> for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: TwistedEdwardsCurve<F>) -> Self {
        TwistedEdwardsMontgomeryReduction::from_twisted_edwards(curve).montgomery_curve
    }
}

impl<F: Field> From<&MontgomeryCurve<F>> for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: &MontgomeryCurve<F>) -> Self {
        TwistedEdwardsMontgomeryReduction::from_montgomery(curve.clone()).twisted_edwards_curve
    }
}

impl<F: Field> From<MontgomeryCurve<F>> for TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: MontgomeryCurve<F>) -> Self {
        TwistedEdwardsMontgomeryReduction::from_montgomery(curve).twisted_edwards_curve
    }
}
