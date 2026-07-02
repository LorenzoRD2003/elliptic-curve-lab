use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve, MontgomeryCurve, ShortWeierstrassCurve,
    traits::{CurveModel, CurveModelConversion, CurveModelConversionError},
};
use crate::fields::traits::*;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};

/// Explicit reduction data from a Montgomery model to a short-Weierstrass
/// companion in characteristic different from `2` and `3`.
///
/// The stored coordinate change follows the affine formulas
///
/// `x = B X - A / 3`, and `y = B Y`
///
/// from short coordinates `(X, Y)` back to Montgomery coordinates `(x, y)`.
pub(crate) struct MontgomeryShortReduction<F: Field> {
    montgomery_curve: MontgomeryCurve<F>,
    short_curve: ShortWeierstrassCurve<F>,
    shift: F::Elem,
}

impl<F: Field> MontgomeryShortReduction<F> {
    /// Builds the explicit reduction object attached to one Montgomery curve.
    pub(crate) fn new(montgomery_curve: MontgomeryCurve<F>) -> Result<Self, CurveError> {
        if F::has_characteristic(2) || F::has_characteristic(3) {
            let characteristic = F::characteristic().to_biguint();
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let one_third = invert_small_integer::<F>(3);
        let shift = F::mul(montgomery_curve.a(), &one_third);

        let short_curve = ShortWeierstrassCurve::new(
            F::div(
                &F::sub(&F::from_i64(3), &F::square(montgomery_curve.a())),
                &F::mul(&F::from_i64(3), &F::square(montgomery_curve.b())),
            )
            .expect("validated Montgomery curve has B != 0"),
            F::div(
                &F::sub(
                    &F::mul(&F::from_i64(2), &F::cube(montgomery_curve.a())),
                    &F::mul(&F::from_i64(9), montgomery_curve.a()),
                ),
                &F::mul(&F::from_i64(27), &F::cube(montgomery_curve.b())),
            )
            .expect("validated Montgomery curve has B != 0"),
        )?;

        Ok(Self {
            montgomery_curve,
            short_curve,
            shift,
        })
    }

    /// Transports one point from the Montgomery model to its short companion.
    fn transport_point_to_short(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.montgomery_curve.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let numerator = F::add(x, &self.shift);
                let short_x = F::div(&numerator, self.montgomery_curve.b())
                    .expect("validated Montgomery curve has B != 0");
                let short_y = F::div(y, self.montgomery_curve.b())
                    .expect("validated Montgomery curve has B != 0");
                let image = AffinePoint::new(short_x, short_y);
                if !self.short_curve.contains(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }

                Ok(image)
            }
        }
    }

    /// Transports one point from the short companion back to the Montgomery
    /// model.
    fn transport_point_to_montgomery(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.short_curve.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let montgomery_x = F::sub(&F::mul(self.montgomery_curve.b(), x), &self.shift);
                let montgomery_y = F::mul(self.montgomery_curve.b(), y);
                let image = AffinePoint::new(montgomery_x, montgomery_y);
                if !self.montgomery_curve.contains_affine_point(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }

                Ok(image)
            }
        }
    }
}

impl<F: Field> CurveModelConversion for MontgomeryShortReduction<F>
where
    F::Elem: Clone,
{
    type Source = MontgomeryCurve<F>;
    type Target = ShortWeierstrassCurve<F>;

    fn source(&self) -> &Self::Source {
        &self.montgomery_curve
    }

    fn target(&self) -> &Self::Target {
        &self.short_curve
    }

    fn map_source_point(
        &self,
        point: &<Self::Source as CurveModel>::Point,
    ) -> Result<<Self::Target as CurveModel>::Point, CurveModelConversionError> {
        self.transport_point_to_short(point)
            .map_err(|error| match error {
                CurveError::PointNotOnCurve => CurveModelConversionError::PointNotOnSource,
                other => CurveModelConversionError::from(other),
            })
    }

    fn map_target_point(
        &self,
        point: &<Self::Target as CurveModel>::Point,
    ) -> Result<<Self::Source as CurveModel>::Point, CurveModelConversionError> {
        self.transport_point_to_montgomery(point)
            .map_err(|error| match error {
                CurveError::PointNotOnCurve => CurveModelConversionError::PointNotOnTarget,
                other => CurveModelConversionError::from(other),
            })
    }
}

fn invert_small_integer<F: Field>(value: i64) -> F::Elem {
    F::inv(&F::from_i64(value)).unwrap_or_else(|| {
        panic!("small integer {value} should be invertible in this characteristic")
    })
}

impl<F: Field> MontgomeryCurve<F> {
    /// Builds the explicit conversion witness to short-Weierstrass form.
    pub fn conversion_to_short_weierstrass(
        &self,
    ) -> Result<
        impl CurveModelConversion<Source = Self, Target = ShortWeierstrassCurve<F>>,
        CurveModelConversionError,
    >
    where
        F::Elem: Clone,
    {
        MontgomeryShortReduction::new(self.clone()).map_err(Into::into)
    }

    /// Returns the reduced short-Weierstrass companion when the characteristic
    /// permits the classical Montgomery reduction.
    pub fn try_as_short_weierstrass(&self) -> Result<ShortWeierstrassCurve<F>, CurveError> {
        ShortWeierstrassCurve::try_from(self).map_err(CurveError::from)
    }

    /// Returns the same curve viewed inside the general Weierstrass family.
    pub fn as_general_weierstrass(&self) -> GeneralWeierstrassCurve<F>
    where
        F::Elem: Clone,
    {
        GeneralWeierstrassCurve::from(self)
    }
}

impl<F: Field> TryFrom<&MontgomeryCurve<F>> for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: &MontgomeryCurve<F>) -> Result<Self, Self::Error> {
        curve
            .conversion_to_short_weierstrass()
            .map(|conversion| conversion.target().clone())
    }
}

impl<F: Field> TryFrom<MontgomeryCurve<F>> for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: MontgomeryCurve<F>) -> Result<Self, Self::Error> {
        MontgomeryShortReduction::new(curve)
            .map_err(Into::into)
            .map(|conversion| conversion.target().clone())
    }
}

impl<F: Field + EnumerableFiniteField + SqrtField> TryFrom<&ShortWeierstrassCurve<F>>
    for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: &ShortWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        for x in F::elements() {
            let candidate_two_torsion = AffinePoint::new(x.clone(), F::zero());
            if !curve.contains(&candidate_two_torsion) {
                continue;
            }

            let tangent_factor = F::add(&F::mul(&F::from_i64(3), &F::square(&x)), curve.a());
            let Some(u_squared_root) = F::sqrt(&tangent_factor) else {
                continue;
            };
            if F::is_zero(&u_squared_root) {
                continue;
            }

            let b = F::inv(&u_squared_root).ok_or(CurveModelConversionError::Curve(
                CurveError::NoMontgomeryModelOverBaseField,
            ))?;
            let a = F::mul(&F::from_i64(3), &F::mul(&x, &b));

            return MontgomeryCurve::new(a, b).map_err(Into::into);
        }

        Err(CurveModelConversionError::Curve(
            CurveError::NoMontgomeryModelOverBaseField,
        ))
    }
}

impl<F: Field + EnumerableFiniteField + SqrtField> TryFrom<ShortWeierstrassCurve<F>>
    for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: ShortWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        MontgomeryCurve::try_from(&curve)
    }
}

impl<F: Field + EnumerableFiniteField + SqrtField> GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    /// Returns a Montgomery companion over the same base field when the
    /// general curve reduces to a short-Weierstrass model that is itself
    /// compatible with the Montgomery family.
    pub fn try_as_montgomery(&self) -> Result<MontgomeryCurve<F>, CurveModelConversionError> {
        MontgomeryCurve::try_from(self)
    }
}

impl<F: Field + EnumerableFiniteField + SqrtField> TryFrom<&GeneralWeierstrassCurve<F>>
    for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: &GeneralWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        let short_curve = ShortWeierstrassCurve::try_from(curve)?;
        MontgomeryCurve::try_from(&short_curve)
    }
}

impl<F: Field + EnumerableFiniteField + SqrtField> TryFrom<GeneralWeierstrassCurve<F>>
    for MontgomeryCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: GeneralWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        MontgomeryCurve::try_from(&curve)
    }
}

impl<F: Field> From<&MontgomeryCurve<F>> for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: &MontgomeryCurve<F>) -> Self {
        GeneralWeierstrassCurve::new(
            F::zero(),
            F::div(curve.a(), curve.b()).expect("validated Montgomery curve has B != 0"),
            F::zero(),
            F::div(&F::one(), &F::square(curve.b()))
                .expect("validated Montgomery curve has B != 0"),
            F::zero(),
        )
        .expect("a validated Montgomery curve should stay non-singular in the general family")
    }
}

impl<F: Field> From<MontgomeryCurve<F>> for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: MontgomeryCurve<F>) -> Self {
        GeneralWeierstrassCurve::from(&curve)
    }
}
