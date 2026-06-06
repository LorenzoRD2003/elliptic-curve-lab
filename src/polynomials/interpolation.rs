use crate::fields::Field;
use crate::polynomials::{DensePolynomial, PolynomialError};

/// Interpolates a univariate polynomial over a field from sample points using
/// the classical Lagrange formula.
///
/// Given distinct sample points
///
/// `(x_0, y_0), (x_1, y_1), ..., (x_{n-1}, y_{n-1})`
///
/// this function constructs the unique polynomial `p(x)` of degree at most
/// `n - 1` such that
///
/// `p(x_i) = y_i`
///
/// for every sample index `i`.
///
/// The implementation is intentionally direct and educational:
///
/// - it builds each Lagrange basis polynomial explicitly
/// - it uses the already available dense polynomial arithmetic
/// - it favors clarity over asymptotic efficiency
///
/// Behavior notes:
///
/// - an empty sample list returns the zero polynomial
/// - a single sample returns the corresponding constant polynomial
/// - repeated `x` coordinates are rejected, because interpolation over a field
///   requires distinct abscissas in this formulation
///
/// TODO:
///
/// - add Newton interpolation for an incremental formulation
/// - add barycentric interpolation for numerically lighter repeated evaluation
/// - add subproduct-tree / divide-and-conquer interpolation once performance
///   becomes relevant
pub fn lagrange_interpolate<F: Field>(
    samples: &[(F::Elem, F::Elem)],
) -> Result<DensePolynomial<F>, PolynomialError> {
    if samples.is_empty() {
        return Ok(DensePolynomial::new(Vec::new()));
    }

    let mut result = DensePolynomial::<F>::new(Vec::new());

    for (i, (x_i, y_i)) in samples.iter().enumerate() {
        let mut numerator = DensePolynomial::<F>::new(vec![F::one()]);
        let mut denominator = F::one();

        for (j, (x_j, _)) in samples.iter().enumerate() {
            if i == j {
                continue;
            }

            if F::eq(x_i, x_j) {
                return Err(PolynomialError::DuplicateInterpolationAbscissa);
            }

            let factor = DensePolynomial::<F>::new(vec![F::neg(x_j), F::one()]);
            numerator = numerator.mul(&factor);

            let difference = F::sub(x_i, x_j);
            denominator = F::mul(&denominator, &difference);
        }

        let scaling = F::div(y_i, &denominator)
            .map_err(|_| PolynomialError::NonInvertibleInterpolationDenominator)?;
        let scaled_basis = numerator.scale(&scaling);
        result = result.add(&scaled_basis);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use crate::fields::{Field, Fp, Q};
    use crate::polynomials::{DensePolynomial, PolynomialError, evaluation::evaluate_dense};
    use crate::proptest_support::{dense_polynomial, distinct_fp_elements};

    use crate::polynomials::interpolation::lagrange_interpolate;

    type F17 = Fp<17>;
    type F17Samples = Vec<(<F17 as Field>::Elem, <F17 as Field>::Elem)>;

    fn q(numerator: i64, denominator: i64) -> <Q as Field>::Elem {
        let numerator = Q::from_i64(numerator);
        let denominator = Q::from_i64(denominator);
        Q::div(&numerator, &denominator).expect("denominator should be non-zero")
    }

    fn assert_dense_eq<F: Field>(actual: &DensePolynomial<F>, expected: &DensePolynomial<F>) {
        assert_eq!(actual.coefficients().len(), expected.coefficients().len());

        for (actual, expected) in actual.coefficients().iter().zip(expected.coefficients()) {
            assert!(F::eq(actual, expected));
        }
    }

    #[test]
    fn lagrange_interpolate_returns_zero_for_empty_input() {
        let polynomial = lagrange_interpolate::<F17>(&[]).expect("empty interpolation should work");

        assert!(polynomial.is_zero());
        assert_eq!(polynomial.coefficients(), &[]);
    }

    #[test]
    fn lagrange_interpolate_returns_constant_for_single_sample() {
        let polynomial =
            lagrange_interpolate::<F17>(&[(F17::elem_from_u64(9), F17::elem_from_u64(4))])
                .expect("single sample should interpolate");

        assert_dense_eq(
            &polynomial,
            &DensePolynomial::<F17>::new(vec![F17::elem_from_u64(4)]),
        );
    }

    #[test]
    fn lagrange_interpolate_reconstructs_linear_polynomial_over_f17() {
        let samples = [
            (F17::elem_from_u64(0), F17::elem_from_u64(3)),
            (F17::elem_from_u64(1), F17::elem_from_u64(8)),
        ];

        let polynomial = lagrange_interpolate::<F17>(&samples).expect("interpolation should work");

        assert_dense_eq(
            &polynomial,
            &DensePolynomial::<F17>::new(vec![F17::elem_from_u64(3), F17::elem_from_u64(5)]),
        );
    }

    #[test]
    fn lagrange_interpolate_reconstructs_quadratic_polynomial_over_f17() {
        let samples = [
            (F17::elem_from_u64(0), F17::elem_from_u64(3)),
            (F17::elem_from_u64(1), F17::elem_from_u64(10)),
            (F17::elem_from_u64(2), F17::elem_from_u64(4)),
        ];

        let polynomial = lagrange_interpolate::<F17>(&samples).expect("interpolation should work");

        assert_dense_eq(
            &polynomial,
            &DensePolynomial::<F17>::new(vec![
                F17::elem_from_u64(3),
                F17::elem_from_u64(5),
                F17::elem_from_u64(2),
            ]),
        );
    }

    #[test]
    fn lagrange_interpolate_matches_all_input_samples_over_q() {
        let samples = [(q(0, 1), q(1, 2)), (q(1, 1), q(7, 6)), (q(2, 1), q(17, 6))];

        let polynomial = lagrange_interpolate::<Q>(&samples).expect("interpolation should work");

        assert_dense_eq(
            &polynomial,
            &DensePolynomial::<Q>::new(vec![q(1, 2), q(1, 6), q(1, 2)]),
        );

        for (x, y) in &samples {
            let value = evaluate_dense(&polynomial, x).expect("evaluation should work");
            assert!(Q::eq(&value, y));
        }
    }

    #[test]
    fn lagrange_interpolate_rejects_duplicate_x_coordinates() {
        let samples = [
            (F17::elem_from_u64(3), F17::elem_from_u64(1)),
            (F17::elem_from_u64(3), F17::elem_from_u64(9)),
        ];

        let error =
            lagrange_interpolate::<F17>(&samples).expect_err("duplicate x values should fail");

        assert_eq!(error, PolynomialError::DuplicateInterpolationAbscissa);
    }

    fn interpolation_case() -> impl Strategy<Value = (DensePolynomial<F17>, F17Samples)> {
        dense_polynomial::<17>(4).prop_flat_map(|polynomial| {
            let sample_count = polynomial.degree().map_or(1, |degree| degree + 1);
            distinct_fp_elements::<17>(sample_count).prop_map(move |xs| {
                let samples = xs
                    .into_iter()
                    .map(|x| {
                        let y = evaluate_dense(&polynomial, &x).expect("evaluation should succeed");
                        (x, y)
                    })
                    .collect::<Vec<_>>();
                (polynomial.clone(), samples)
            })
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]

        #[test]
        fn property_lagrange_interpolation_recovers_small_dense_polynomials(
            case in interpolation_case(),
        ) {
            let (polynomial, samples) = case;
            let interpolated = lagrange_interpolate::<F17>(&samples).expect("interpolation should succeed");
            prop_assert_eq!(interpolated, polynomial);
        }
    }
}
