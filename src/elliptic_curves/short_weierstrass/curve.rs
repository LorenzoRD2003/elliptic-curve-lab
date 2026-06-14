use core::fmt;

use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::error::CurveError;
use crate::elliptic_curves::isomorphisms::CurveIsomorphismError;
use crate::fields::Field;
use crate::polynomials::DensePolynomial;

/// Short-Weierstrass curve model `y^2 = x^3 + ax + b`.
///
/// This educational implementation currently supports only fields of
/// characteristic different from `2` and `3`, where the classical short form
/// and its discriminant formula behave as expected.
pub struct ShortWeierstrassCurve<F: Field> {
    pub(crate) a: F::Elem,
    pub(crate) b: F::Elem,
}

impl<F: Field> Clone for ShortWeierstrassCurve<F>
where
    F::Elem: Clone,
{
    fn clone(&self) -> Self {
        Self {
            a: self.a.clone(),
            b: self.b.clone(),
        }
    }
}

impl<F: Field> PartialEq for ShortWeierstrassCurve<F>
where
    F::Elem: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b
    }
}

impl<F: Field> Eq for ShortWeierstrassCurve<F> where F::Elem: Eq {}

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Builds a validated short-Weierstrass curve descriptor.
    pub fn new(a: F::Elem, b: F::Elem) -> Result<Self, CurveError> {
        let characteristic = F::characteristic();
        if matches!(characteristic, 2 | 3) {
            return Err(CurveError::UnsupportedCharacteristic { characteristic });
        }

        let curve = Self { a, b };
        if F::is_zero(&curve.discriminant()) {
            return Err(CurveError::SingularCurve);
        }

        Ok(curve)
    }

    /// Returns the `a` coefficient in `x^3 + ax + b`.
    pub fn a(&self) -> &F::Elem {
        &self.a
    }

    /// Returns the `b` coefficient in `x^3 + ax + b`.
    pub fn b(&self) -> &F::Elem {
        &self.b
    }

    /// Returns the short-Weierstrass equation as plain text.
    pub fn to_equation_string(&self) -> String
    where
        F::Elem: fmt::Display,
    {
        format!("y^2 = x^3 + ({})x + ({})", self.a, self.b)
    }

    /// Returns the cubic polynomial `x^3 + ax + b`.
    ///
    /// The coefficients are returned in ascending degree order, so the dense
    /// representation is `[b, a, 0, 1]`.
    pub fn to_cubic(&self) -> DensePolynomial<F> {
        DensePolynomial::<F>::new(vec![
            self.b.clone(),
            self.a.clone(),
            F::zero(),
            F::one(),
        ])
    }

    /// Returns the discriminant `Δ = -16(4a^3 + 27b^2)`.
    pub fn discriminant(&self) -> F::Elem {
        let four = F::from_i64(4);
        let minus_sixteen = F::from_i64(-16);
        let twenty_seven = F::from_i64(27);

        let four_a_cubed = F::mul(&four, &F::cube(&self.a));
        let twenty_seven_b_squared = F::mul(&twenty_seven, &F::square(&self.b));
        let inner = F::add(&four_a_cubed, &twenty_seven_b_squared);
        F::mul(&minus_sixteen, &inner)
    }

    /// Returns the classical Weierstrass invariant `c4 = -48a`.
    pub fn c4(&self) -> F::Elem {
        F::mul(&F::from_i64(-48), &self.a)
    }

    /// Returns the classical Weierstrass invariant `c6 = -864b`.
    pub fn c6(&self) -> F::Elem {
        F::mul(&F::from_i64(-864), &self.b)
    }

    /// Returns the `j`-invariant `j = c4^3 / Δ`.
    pub fn j_invariant(&self) -> F::Elem {
        let c4_cubed = F::cube(&self.c4());
        F::div(&c4_cubed, &self.discriminant())
            .expect("validated short Weierstrass curve has non-zero discriminant")
    }

    /// Returns whether this curve and `other` have the same `j`-invariant.
    pub fn has_same_j_invariant(&self, other: &Self) -> bool {
        F::eq(&self.j_invariant(), &other.j_invariant())
    }

    /// Returns the short-Weierstrass model obtained from the scaling
    /// parameter `u`.
    pub fn scaled_by(&self, u: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&u).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }

        let u2 = F::square(&u);
        let u4 = F::square(&u2);
        let u6 = F::mul(&u4, &u2);

        Self::new(F::mul(&u4, &self.a), F::mul(&u6, &self.b)).map_err(Into::into)
    }

    /// Returns whether `other` is exactly the short-Weierstrass model obtained
    /// by scaling this curve with the supplied parameter `u`.
    pub fn isomorphic_via_scale(&self, other: &Self, u: &F::Elem) -> bool {
        match self.scaled_by(u.clone()) {
            Ok(scaled_curve) => {
                F::eq(scaled_curve.a(), other.a()) && F::eq(scaled_curve.b(), other.b())
            }
            Err(_) => false,
        }
    }

    /// Returns the quadratic twist determined by the non-zero factor `d`.
    pub fn quadratic_twist(&self, d: F::Elem) -> Result<Self, CurveIsomorphismError> {
        if F::inv(&d).is_none() {
            return Err(CurveIsomorphismError::NonInvertibleScale);
        }

        let d2 = F::square(&d);
        let d3 = F::mul(&d2, &d);

        Self::new(F::mul(&d2, &self.a), F::mul(&d3, &self.b)).map_err(Into::into)
    }

    /// Builds a finite affine point without checking the curve equation.
    pub(crate) fn unchecked_point(&self, x: F::Elem, y: F::Elem) -> AffinePoint<F> {
        AffinePoint::new(x, y)
    }

    /// Returns the cubic right-hand side `x^3 + ax + b`.
    pub(crate) fn rhs_value(&self, x: &F::Elem) -> F::Elem {
        let x_cubed = F::cube(x);
        let ax = F::mul(&self.a, x);
        F::add(&F::add(&x_cubed, &ax), &self.b)
    }
}
