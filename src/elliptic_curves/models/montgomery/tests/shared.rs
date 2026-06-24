use num_bigint::BigInt;
use num_rational::BigRational;

use crate::elliptic_curves::{AffinePoint, CurveError, MontgomeryCurve, NormalizedMontgomeryCurve};
use crate::fields::{
    Fp,
    traits::{Field, SqrtField},
};

pub(super) type F2 = Fp<2>;
pub(super) type F3 = Fp<3>;
pub(super) type F5 = Fp<5>;
pub(super) type F7 = Fp<7>;

pub(super) fn q(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

pub(super) fn f3_curve() -> MontgomeryCurve<F3> {
    MontgomeryCurve::<F3>::new(F3::zero(), F3::one()).expect("valid Montgomery curve")
}

pub(super) fn f5_curve() -> MontgomeryCurve<F5> {
    MontgomeryCurve::<F5>::new(F5::one(), F5::one()).expect("valid Montgomery curve")
}

pub(super) fn f7_curve() -> MontgomeryCurve<F7> {
    MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::one()).expect("valid Montgomery curve")
}

pub(super) fn f7_scaled_curve() -> MontgomeryCurve<F7> {
    MontgomeryCurve::<F7>::new(F7::from_i64(3), F7::from_i64(2)).expect("valid Montgomery curve")
}

pub(super) fn f7_nonsquare_scaled_curve() -> MontgomeryCurve<F7> {
    MontgomeryCurve::<F7>::new(F7::one(), F7::from_i64(3)).expect("valid Montgomery curve")
}

pub(super) fn normalization_sqrt_b<F>(curve: &MontgomeryCurve<F>) -> F::Elem
where
    F: Field + SqrtField,
{
    F::sqrt(curve.b()).expect("test helper requires B to be a square over the base field")
}

pub(super) fn normalize_point<F>(
    curve: &MontgomeryCurve<F>,
    normalized: &NormalizedMontgomeryCurve<F>,
    point: &AffinePoint<F>,
) -> Result<AffinePoint<F>, CurveError>
where
    F: Field + SqrtField,
    F::Elem: Clone,
{
    if !curve.contains_affine_point(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    match point {
        AffinePoint::Infinity => Ok(AffinePoint::Infinity),
        AffinePoint::Finite { x, y } => {
            let image = AffinePoint::new(x.clone(), F::mul(&normalization_sqrt_b(curve), y));
            if !normalized.contains_affine_point(&image) {
                return Err(CurveError::PointNotOnCurve);
            }
            Ok(image)
        }
    }
}

pub(super) fn denormalize_point<F>(
    curve: &MontgomeryCurve<F>,
    normalized: &NormalizedMontgomeryCurve<F>,
    point: &AffinePoint<F>,
) -> Result<AffinePoint<F>, CurveError>
where
    F: Field + SqrtField,
    F::Elem: Clone,
{
    if !normalized.contains_affine_point(point) {
        return Err(CurveError::PointNotOnCurve);
    }

    match point {
        AffinePoint::Infinity => Ok(AffinePoint::Infinity),
        AffinePoint::Finite { x, y } => {
            let image = AffinePoint::new(
                x.clone(),
                F::div(y, &normalization_sqrt_b(curve))
                    .expect("sqrt(B) is nonzero on a validated Montgomery curve"),
            );
            if !curve.contains_affine_point(&image) {
                return Err(CurveError::PointNotOnCurve);
            }
            Ok(image)
        }
    }
}
