use crate::elliptic_curves::short_weierstrass::schoof::{
    ReducedCurveFunction, ReducedCurveQuotient,
};
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

/// One reduced map of the form `(a(x), b(x) y)` over
/// `F[x,y] / (y^2 - f(x), g(x))`.
///
/// This type models the restricted coordinate shape that appears in Schoof's
/// odd-prime torsion arithmetic. It should be read as the reduced map
///
/// - `x ↦ a(x)`
/// - `y ↦ b(x) y`
///
/// where both `a(x)` and `b(x)` are stored modulo `g(x)`.
#[derive(Debug)]
pub struct ReducedEndomorphism<F: FiniteField> {
    x_map: DensePolynomial<F>,
    y_scale: DensePolynomial<F>,
}

impl<F: FiniteField> Clone for ReducedEndomorphism<F> {
    fn clone(&self) -> Self {
        Self {
            x_map: self.x_map.clone(),
            y_scale: self.y_scale.clone(),
        }
    }
}

impl<F: FiniteField> PartialEq for ReducedEndomorphism<F> {
    fn eq(&self, other: &Self) -> bool {
        self.x_map == other.x_map && self.y_scale == other.y_scale
    }
}

impl<F: FiniteField> ReducedEndomorphism<F> {
    /// Builds one reduced map `(a(x), b(x) y)`.
    ///
    /// Both stored polynomials are reduced modulo the active univariate
    /// modulus `g(x)` before storage.
    ///
    /// Complexity: if `m = deg g`, `a = deg x_map`, and `b = deg y_scale`,
    /// then under the current dense backend this constructor costs
    /// `Θ(m(a + b))` field operations.
    pub(crate) fn new(
        quotient: &ReducedCurveQuotient<F>,
        x_map: DensePolynomial<F>,
        y_scale: DensePolynomial<F>,
    ) -> Self {
        Self {
            x_map: quotient.reduce_poly(&x_map),
            y_scale: quotient.reduce_poly(&y_scale),
        }
    }

    /// Returns the identity reduced map `x ↦ x`, `y ↦ y`.
    ///
    /// Complexity: `Θ(deg g)` field operations.
    pub(crate) fn identity(quotient: &ReducedCurveQuotient<F>) -> Self {
        Self::new(
            quotient,
            DensePolynomial::new(vec![F::zero(), F::one()]),
            DensePolynomial::constant(F::one()),
        )
    }

    /// Returns the stored reduced polynomial `a(x)` in the map
    /// `(a(x), b(x) y)`.
    pub fn x_map(&self) -> &DensePolynomial<F> {
        &self.x_map
    }

    /// Returns the stored reduced polynomial `b(x)` in the map
    /// `(a(x), b(x) y)`.
    pub fn y_scale(&self) -> &DensePolynomial<F> {
        &self.y_scale
    }

    /// Returns the additive inverse in the reduced affine chart.
    ///
    /// This keeps the same `x`-map and negates only the `y`-scale:
    ///
    /// `(a(x), b(x) y) ↦ (a(x), -b(x) y)`.
    ///
    /// Because `y_scale` is already reduced modulo `g(x)`, its negation is
    /// too, so no extra reduction is needed.
    ///
    /// Complexity: `Θ(m)` field operations on a degree-`< m` representative.
    pub(crate) fn additive_inverse(&self) -> Self {
        Self {
            x_map: self.x_map.clone(),
            y_scale: self.y_scale.neg(),
        }
    }

    /// Applies the reduced map to the distinguished `x` coordinate.
    ///
    /// The resulting quotient value is the class `a(x) + y * 0`.
    ///
    /// Complexity: if `m = deg g` and `a = deg x_map`, then under the current
    /// dense backend this costs `Θ(ma)` field operations.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn apply_to_x(&self, quotient: &ReducedCurveQuotient<F>) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            quotient,
            self.x_map.clone(),
            DensePolynomial::new(Vec::new()),
        )
    }

    /// Applies the reduced map to the distinguished `y` coordinate.
    ///
    /// The resulting quotient value is the class `0 + y b(x)`.
    ///
    /// Complexity: if `m = deg g` and `b = deg y_scale`, then under the
    /// current dense backend this costs `Θ(mb)` field operations.
    #[cfg_attr(not(test), allow(dead_code))]
    pub(crate) fn apply_to_y(&self, quotient: &ReducedCurveQuotient<F>) -> ReducedCurveFunction<F> {
        ReducedCurveFunction::new(
            quotient,
            DensePolynomial::new(Vec::new()),
            self.y_scale.clone(),
        )
    }

    /// Evaluates one outer polynomial at the stored reduced `x`-image and
    /// reduces the result modulo `g(x)`.
    ///
    /// This uses Horner's rule directly in the quotient `F[x] / (g(x))`.
    /// Concretely, if
    ///
    /// `outer(x) = c_0 + c_1 x + ... + c_d x^d`
    ///
    /// and the stored reduced `x`-image is `a(x)`, then this method computes
    /// the class of `outer(a(x))` through the recurrence
    ///
    /// `value <- value * a(x) + c_i`
    ///
    /// iterating from the top coefficient downward, with a reduction modulo
    /// `g(x)` after each multiplication and addition step.
    ///
    /// This is preferable to first expanding `outer(a(x))` as one large
    /// polynomial in `F[x]` and reducing only at the end: the Horner route
    /// keeps every intermediate representative inside the quotient
    /// `F[x] / (g(x))`, so degrees stay controlled by `deg g` rather than
    /// growing toward `deg(outer) * deg(a)`.
    ///
    /// Complexity: if `d = deg outer` and `m = deg g`, then under the current
    /// dense backend this costs `Θ(d m^2)` field operations.
    pub(crate) fn evaluate_x_at_map_mod(
        &self,
        quotient: &ReducedCurveQuotient<F>,
        outer: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let mut value = DensePolynomial::new(Vec::new());

        for coefficient in outer.coefficients().iter().rev() {
            value = quotient.reduce_poly(&value.mul(&self.x_map));
            value =
                quotient.reduce_poly(&value.add(&DensePolynomial::constant(coefficient.clone())));
        }
        value
    }

    /// Composes two reduced maps of the form `(a(x), b(x) y)`.
    ///
    /// In the current Schoof arithmetic, this composition is the
    /// multiplicative operation on the represented reduced endomorphisms.
    ///
    /// If `self = (a1(x), b1(x) y)` and `rhs = (a2(x), b2(x) y)`, the
    /// composition `self ∘ rhs` is
    ///
    /// - `x ↦ a1(a2(x))`
    /// - `y ↦ b1(a2(x)) b2(x) y`
    ///
    /// with all polynomial data reduced modulo `g(x)`.
    ///
    /// Complexity: under the current dense backend, composition is dominated
    /// by two substitutions of degree-`< m` polynomials into a degree-`< m`
    /// representative, so it costs `Θ(m^3)` field operations when `m = deg g`.
    pub(crate) fn compose(&self, quotient: &ReducedCurveQuotient<F>, rhs: &Self) -> Self {
        let x_map = rhs.evaluate_x_at_map_mod(quotient, &self.x_map);
        let y_scale = quotient.reduce_poly(
            &rhs.evaluate_x_at_map_mod(quotient, &self.y_scale)
                .mul(&rhs.y_scale),
        );

        Self::new(quotient, x_map, y_scale)
    }
}

#[cfg(test)]
mod tests {
    use crate::elliptic_curves::{
        ShortWeierstrassCurve,
        short_weierstrass::schoof::{ReducedCurveQuotient, ReducedEndomorphism},
    };
    use crate::fields::{Fp, traits::Field};
    use crate::polynomials::DensePolynomial;

    type F7 = Fp<7>;

    fn sample_curve() -> ShortWeierstrassCurve<F7> {
        ShortWeierstrassCurve::<F7>::new(F7::from_i64(2), F7::from_i64(3))
            .expect("sample curve should be smooth")
    }

    fn sample_quotient() -> ReducedCurveQuotient<F7> {
        ReducedCurveQuotient::new(
            sample_curve(),
            DensePolynomial::new(vec![F7::one(), F7::zero(), F7::one()]),
        )
        .expect("non-zero modulus should define a quotient")
    }

    #[test]
    fn identity_endomorphism_has_the_expected_coordinate_maps() {
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);

        assert_eq!(
            identity.x_map(),
            &DensePolynomial::new(vec![F7::zero(), F7::one()])
        );
        assert_eq!(identity.y_scale(), &DensePolynomial::constant(F7::one()));
    }

    #[test]
    fn apply_to_x_and_y_return_the_expected_quotient_values() {
        let quotient = sample_quotient();
        let endomorphism = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::from_i64(5), F7::from_i64(1)]),
        );

        let x_image = endomorphism.apply_to_x(&quotient);
        let y_image = endomorphism.apply_to_y(&quotient);

        assert_eq!(x_image.x_part(), endomorphism.x_map());
        assert!(x_image.y_part().is_zero());
        assert!(y_image.x_part().is_zero());
        assert_eq!(y_image.y_part(), endomorphism.y_scale());
    }

    #[test]
    fn substitution_uses_the_stored_x_image_modulo_g() {
        let quotient = sample_quotient();
        let endomorphism = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(1)]),
            DensePolynomial::constant(F7::one()),
        );
        let outer = DensePolynomial::new(vec![F7::from_i64(2), F7::from_i64(3), F7::from_i64(4)]);

        let substituted = endomorphism.evaluate_x_at_map_mod(&quotient, &outer);
        // Here `g(x) = x^2 + 1`, so in the quotient we have `x^2 = -1`.
        // For `a(x) = 1 + x` and `outer(x) = 2 + 3x + 4x^2`,
        // `outer(a(x)) = 2 + 3(1 + x) + 4(1 + x)^2 = 5 + 4x` modulo `g(x)`.
        let expected = DensePolynomial::new(vec![F7::from_i64(5), F7::from_i64(4)]);

        assert_eq!(substituted, expected);
    }

    #[test]
    fn composition_matches_the_coordinate_formula() {
        let quotient = sample_quotient();
        let first = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(2), F7::from_i64(1)]),
            DensePolynomial::new(vec![F7::from_i64(4), F7::from_i64(3)]),
        );
        let second = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(6), F7::from_i64(1)]),
            DensePolynomial::new(vec![F7::from_i64(5)]),
        );

        let composed = first.compose(&quotient, &second);
        let expected_x = second.evaluate_x_at_map_mod(&quotient, first.x_map());
        let expected_y = quotient.reduce_poly(
            &second
                .evaluate_x_at_map_mod(&quotient, first.y_scale())
                .mul(second.y_scale()),
        );

        assert_eq!(composed.x_map(), &expected_x);
        assert_eq!(composed.y_scale(), &expected_y);
    }

    #[test]
    fn identity_composes_as_expected_on_both_sides() {
        let quotient = sample_quotient();
        let identity = ReducedEndomorphism::identity(&quotient);
        let endomorphism = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(6)]),
        );

        let left = identity.compose(&quotient, &endomorphism);
        let right = endomorphism.compose(&quotient, &identity);

        assert_eq!(left.x_map(), endomorphism.x_map());
        assert_eq!(left.y_scale(), endomorphism.y_scale());
        assert_eq!(right.x_map(), endomorphism.x_map());
        assert_eq!(right.y_scale(), endomorphism.y_scale());
    }

    #[test]
    fn additive_inverse_keeps_x_map_and_negates_y_scale() {
        let quotient = sample_quotient();
        let endomorphism = ReducedEndomorphism::new(
            &quotient,
            DensePolynomial::new(vec![F7::from_i64(3), F7::from_i64(2)]),
            DensePolynomial::new(vec![F7::from_i64(1), F7::from_i64(6)]),
        );

        let inverse = endomorphism.additive_inverse();

        assert_eq!(inverse.x_map(), endomorphism.x_map());
        assert_eq!(inverse.y_scale(), &endomorphism.y_scale().neg());
    }
}
