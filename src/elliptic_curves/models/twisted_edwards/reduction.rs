use crate::elliptic_curves::{MontgomeryCurve, TwistedEdwardsCurve};
use crate::fields::traits::Field;

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
}

impl<F: Field> TwistedEdwardsCurve<F>
where
    F::Elem: Clone,
{
    /// Returns the canonically associated Montgomery companion
    ///
    /// `E_{A,B}: B y^2 = x^3 + A x^2 + x`
    ///
    /// with `A = 2(a + d) / (a - d)`, `B = 4 / (a - d)`.
    pub fn as_montgomery(&self) -> MontgomeryCurve<F> {
        MontgomeryCurve::from(self)
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
    pub fn as_twisted_edwards(&self) -> TwistedEdwardsCurve<F> {
        TwistedEdwardsCurve::from(self)
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
