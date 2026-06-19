use crate::elliptic_curves::short_weierstrass::schoof::ReducedCurveQuotient;
use crate::fields::traits::FiniteField;
use crate::polynomials::DensePolynomial;

/// Result of attempting to invert one univariate polynomial modulo the stored
/// quotient modulus `g(x)`.
///
/// In the Schoof workflow, failure of invertibility is not an exceptional
/// condition. A non-unit denominator reveals a non-trivial common factor with
/// the active modulus, which is precisely the information the algorithm wants
/// to surface.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum QuotientInverseResult<F: FiniteField> {
    /// The input defines a unit class in `F[x] / (g(x))`, together with one
    /// reduced inverse representative.
    Inverse(DensePolynomial<F>),
    /// The input is not a unit modulo `g(x)`.
    ///
    /// The returned witness is the monic gcd of the input with `g(x)`. In
    /// Schoof's algorithm this gcd is useful data, not just an error report.
    NonUnit { witness_gcd: DensePolynomial<F> },
}

struct ExtendedEuclideanInverseData<F: FiniteField> {
    gcd: DensePolynomial<F>,
    bezout_coefficient_for_input: DensePolynomial<F>,
}

impl<F: FiniteField> ReducedCurveQuotient<F> {
    /// Attempts to invert one univariate polynomial modulo `g(x)`.
    ///
    /// If `poly` is a unit in `F[x] / (g(x))`, this returns one reduced
    /// inverse representative. Otherwise it returns the monic gcd
    /// `gcd(poly, g)`, which witnesses exactly why the class is not a unit.
    ///
    /// The algorithm is the extended Euclidean algorithm in `F[x]`. It does not
    /// return a fallible `Result`, because non-invertibility is expected and
    /// useful in Schoof's quotient arithmetic.
    ///
    /// Complexity: one polynomial reduction plus one extended Euclidean run on
    /// degrees bounded by `m = deg g`. Under the current dense backend this is
    /// `Θ(m^2)` field operations.
    pub(crate) fn try_invert_poly(&self, poly: &DensePolynomial<F>) -> QuotientInverseResult<F> {
        let reduced_poly = self.reduce_poly(poly);
        self.reject_zero_as_non_unit(&reduced_poly)
            .unwrap_or_else(|| {
                self.finish_inverse_attempt(self.extended_euclidean_inverse_data(reduced_poly))
            })
    }

    fn reject_zero_as_non_unit(
        &self,
        reduced_poly: &DensePolynomial<F>,
    ) -> Option<QuotientInverseResult<F>> {
        reduced_poly
            .is_zero()
            .then(|| QuotientInverseResult::NonUnit {
                witness_gcd: self.modulus().clone(),
            })
    }

    fn extended_euclidean_inverse_data(
        &self,
        reduced_poly: DensePolynomial<F>,
    ) -> ExtendedEuclideanInverseData<F> {
        let mut previous_remainder = self.modulus().clone();
        let mut remainder = reduced_poly;
        let mut previous_coefficient = DensePolynomial::new(Vec::new());
        let mut coefficient = DensePolynomial::constant(F::one());

        while !remainder.is_zero() {
            let (quotient, next_remainder) = previous_remainder
                .div_rem(&remainder)
                .expect("euclidean step divides by a non-zero polynomial");
            let next_coefficient = previous_coefficient.sub(&quotient.mul(&coefficient));

            previous_remainder = remainder;
            remainder = next_remainder;
            previous_coefficient = coefficient;
            coefficient = next_coefficient;
        }

        ExtendedEuclideanInverseData {
            gcd: previous_remainder,
            bezout_coefficient_for_input: previous_coefficient,
        }
    }

    fn finish_inverse_attempt(
        &self,
        data: ExtendedEuclideanInverseData<F>,
    ) -> QuotientInverseResult<F> {
        if self.gcd_is_unit(&data.gcd) {
            QuotientInverseResult::Inverse(
                self.normalize_unit_inverse(data.gcd, data.bezout_coefficient_for_input),
            )
        } else {
            QuotientInverseResult::NonUnit {
                witness_gcd: self.normalize_non_unit_gcd(data.gcd),
            }
        }
    }

    fn gcd_is_unit(&self, gcd: &DensePolynomial<F>) -> bool {
        gcd.degree() == Some(0)
    }

    fn normalize_unit_inverse(
        &self,
        gcd: DensePolynomial<F>,
        bezout_coefficient_for_input: DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let scale = F::inverse(
            gcd.leading_coefficient()
                .expect("non-zero constant gcd has a leading coefficient"),
        )
        .expect("non-zero field element should be invertible");
        self.reduce_poly(&bezout_coefficient_for_input.scale(&scale))
    }

    fn normalize_non_unit_gcd(&self, gcd: DensePolynomial<F>) -> DensePolynomial<F> {
        gcd.make_monic()
            .expect("non-zero gcd admits monic normalization")
    }
}
