use crate::elliptic_curves::{
    AffinePoint, CurveError, GeneralWeierstrassCurve, ShortWeierstrassCurve,
    traits::{
        CurveModel, CurveModelConversion, CurveModelConversionError, ReversedCurveModelConversion,
    },
};
use crate::fields::traits::*;

/// Explicit reduction data from a general Weierstrass model to a
/// short-Weierstrass companion in characteristic different from `2` and `3`.
///
/// The stored coordinate change follows the admissible transformation
///
/// `x = X + r`, and `y = Y + sX + t`
///
/// from short coordinates `(X, Y)` back to the original general coordinates
/// `(x, y)`.
pub(crate) struct GeneralWeierstrassReduction<F: Field> {
    general_curve: GeneralWeierstrassCurve<F>,
    short_curve: ShortWeierstrassCurve<F>,
    r: F::Elem,
    s: F::Elem,
    t: F::Elem,
}

impl<F: Field> GeneralWeierstrassReduction<F> {
    /// Builds the explicit reduction object attached to one general
    /// Weierstrass curve.
    pub(crate) fn new(general_curve: GeneralWeierstrassCurve<F>) -> Result<Self, CurveError> {
        if F::has_characteristic(2) || F::has_characteristic(3) {
            let characteristic = F::characteristic().to_biguint();
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let r = F::neg(&F::mul(&general_curve.b2(), &invert_small_integer::<F>(12)));
        let s = F::neg(&F::mul(general_curve.a1(), &invert_small_integer::<F>(2)));
        let t = F::neg(&F::mul(
            &F::add(general_curve.a3(), &F::mul(&r, general_curve.a1())),
            &invert_small_integer::<F>(2),
        ));

        let short_curve = ShortWeierstrassCurve::new(
            F::neg(&F::mul(&general_curve.c4(), &invert_small_integer::<F>(48))),
            F::neg(&F::mul(
                &general_curve.c6(),
                &invert_small_integer::<F>(864),
            )),
        )?;

        Ok(Self {
            general_curve,
            short_curve,
            r,
            s,
            t,
        })
    }

    /// Transports one point from the general model to its short companion.
    fn transport_point_to_short(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.general_curve.contains_affine_point(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let short_x = F::sub(x, &self.r);
                let short_y = F::sub(y, &F::add(&F::mul(&self.s, &short_x), &self.t));
                let image = AffinePoint::new(short_x, short_y);
                if !self.short_curve.contains(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }

                Ok(image)
            }
        }
    }

    /// Transports one point from the short companion back to the general
    /// Weierstrass model.
    fn transport_point_to_general(
        &self,
        point: &AffinePoint<F>,
    ) -> Result<AffinePoint<F>, CurveError> {
        if !self.short_curve.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }

        match point {
            AffinePoint::Infinity => Ok(AffinePoint::Infinity),
            AffinePoint::Finite { x, y } => {
                let general_x = F::add(x, &self.r);
                let general_y = F::add(y, &F::add(&F::mul(&self.s, x), &self.t));
                let image = AffinePoint::new(general_x, general_y);
                if !self.general_curve.contains_affine_point(&image) {
                    return Err(CurveError::PointNotOnCurve);
                }

                Ok(image)
            }
        }
    }
}

impl<F: Field> GeneralWeierstrassReduction<F>
where
    F::Elem: Clone,
{
    /// Builds the explicit inclusion of a short-Weierstrass curve into the
    /// general Weierstrass family.
    pub(crate) fn from_short_weierstrass(short_curve: ShortWeierstrassCurve<F>) -> Self {
        let general_curve = GeneralWeierstrassCurve::new(
            F::zero(),
            F::zero(),
            F::zero(),
            short_curve.a().clone(),
            short_curve.b().clone(),
        )
        .expect(
            "a validated short-Weierstrass curve should stay non-singular in the general family",
        );

        Self {
            general_curve,
            short_curve,
            r: F::zero(),
            s: F::zero(),
            t: F::zero(),
        }
    }
}

impl<F: Field> CurveModelConversion for GeneralWeierstrassReduction<F>
where
    F::Elem: Clone,
{
    type Source = GeneralWeierstrassCurve<F>;
    type Target = ShortWeierstrassCurve<F>;

    fn source(&self) -> &Self::Source {
        &self.general_curve
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
        self.transport_point_to_general(point)
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

impl<F: Field> GeneralWeierstrassCurve<F> {
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
        GeneralWeierstrassReduction::new(self.clone()).map_err(Into::into)
    }

    /// Returns the reduced short-Weierstrass companion when the characteristic
    /// permits the classical reduction.
    pub fn try_as_short_weierstrass(&self) -> Result<ShortWeierstrassCurve<F>, CurveError> {
        ShortWeierstrassCurve::try_from(self).map_err(CurveError::from)
    }
}

impl<F: Field> ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    /// Builds the explicit inclusion of this short-Weierstrass curve into the
    /// general Weierstrass family.
    pub fn conversion_to_general_weierstrass(
        &self,
    ) -> impl CurveModelConversion<Source = Self, Target = GeneralWeierstrassCurve<F>> {
        ReversedCurveModelConversion(GeneralWeierstrassReduction::from_short_weierstrass(
            self.clone(),
        ))
    }

    /// Returns the same curve viewed inside the general Weierstrass family.
    pub fn as_general_weierstrass(&self) -> GeneralWeierstrassCurve<F> {
        GeneralWeierstrassCurve::from(self)
    }
}

impl<F: Field> TryFrom<&GeneralWeierstrassCurve<F>> for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: &GeneralWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        curve
            .conversion_to_short_weierstrass()
            .map(|conversion| conversion.target().clone())
    }
}

impl<F: Field> TryFrom<GeneralWeierstrassCurve<F>> for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    type Error = CurveModelConversionError;

    fn try_from(curve: GeneralWeierstrassCurve<F>) -> Result<Self, Self::Error> {
        GeneralWeierstrassReduction::new(curve)
            .map_err(Into::into)
            .map(|conversion| conversion.target().clone())
    }
}

impl<F: Field> From<&ShortWeierstrassCurve<F>> for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: &ShortWeierstrassCurve<F>) -> Self {
        curve.conversion_to_general_weierstrass().target().clone()
    }
}

impl<F: Field> From<ShortWeierstrassCurve<F>> for GeneralWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn from(curve: ShortWeierstrassCurve<F>) -> Self {
        GeneralWeierstrassReduction::from_short_weierstrass(curve).general_curve
    }
}
