use crate::polynomials::{
    DensePolynomial, PolynomialError,
    irreducibility::{IrreducibilityBackend, IrreducibilityStatus},
};

type DenseIrreducibilityWitness<F> = (DensePolynomial<F>, DensePolynomial<F>);

impl<F> IrreducibilityBackend for F
where
    F: FiniteField + EnumerableFiniteField,
{
    fn irreducibility_status_impl(
        polynomial: &DensePolynomial<Self>,
    ) -> Result<IrreducibilityStatus<Self>, PolynomialError> {
        Self::check_structure().map_err(PolynomialError::InvalidBaseField)?;

        match polynomial.degree() {
            None | Some(0) => return Ok(IrreducibilityStatus::Constant),
            Some(1) => return Ok(IrreducibilityStatus::Linear),
            Some(_) => {}
        }

        let original_leading = polynomial
            .leading_coefficient()
            .expect("non-constant polynomial has a leading coefficient")
            .clone();
        let monic_polynomial = polynomial.make_monic()?;
        let degree = monic_polynomial
            .degree()
            .expect("non-constant polynomial has a degree");

        for divisor_degree in 1..=(degree / 2) {
            if let Some((divisor, quotient)) =
                find_monic_divisor::<Self>(&monic_polynomial, &original_leading, divisor_degree)?
            {
                return Ok(IrreducibilityStatus::Reducible { divisor, quotient });
            }
        }

        Ok(IrreducibilityStatus::Irreducible)
    }
}

fn find_monic_divisor<F>(
    monic_polynomial: &DensePolynomial<F>,
    original_leading: &F::Elem,
    divisor_degree: usize,
) -> Result<Option<DenseIrreducibilityWitness<F>>, PolynomialError>
where
    F: FiniteField + EnumerableFiniteField,
{
    let mut coefficients = Vec::with_capacity(divisor_degree);
    search_monic_divisor_coefficients::<F>(
        monic_polynomial,
        original_leading,
        divisor_degree,
        0,
        &mut coefficients,
    )
}

fn search_monic_divisor_coefficients<F>(
    monic_polynomial: &DensePolynomial<F>,
    original_leading: &F::Elem,
    divisor_degree: usize,
    next_degree: usize,
    coefficients: &mut Vec<F::Elem>,
) -> Result<Option<DenseIrreducibilityWitness<F>>, PolynomialError>
where
    F: FiniteField + EnumerableFiniteField,
{
    if next_degree == divisor_degree {
        let mut divisor_coefficients = coefficients.clone();
        divisor_coefficients.push(F::one());
        let divisor = DensePolynomial::<F>::new(divisor_coefficients);
        let (quotient_for_monic_polynomial, remainder) = monic_polynomial.div_rem(&divisor)?;

        if remainder.is_zero() {
            let quotient = quotient_for_monic_polynomial.scale(original_leading);
            return Ok(Some((divisor, quotient)));
        }

        return Ok(None);
    }

    for value in F::elements() {
        coefficients.push(value);

        if let Some(witness) = search_monic_divisor_coefficients::<F>(
            monic_polynomial,
            original_leading,
            divisor_degree,
            next_degree + 1,
            coefficients,
        )? {
            return Ok(Some(witness));
        }

        coefficients.pop();
    }

    Ok(None)
}
use crate::fields::traits::{EnumerableFiniteField, FiniteField};
