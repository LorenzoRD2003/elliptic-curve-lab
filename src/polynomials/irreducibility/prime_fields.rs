use crate::fields::{Field, FiniteField, Fp};
use crate::polynomials::{DensePolynomial, PolynomialError};

use crate::polynomials::irreducibility::{IrreducibilityBackend, IrreducibilityStatus};

type DenseIrreducibilityWitness<const P: u64> = (DensePolynomial<Fp<P>>, DensePolynomial<Fp<P>>);

impl<const P: u64> IrreducibilityBackend for Fp<P> {
    fn irreducibility_status_impl(
        polynomial: &DensePolynomial<Self>,
    ) -> Result<IrreducibilityStatus<Self>, PolynomialError> {
        Fp::<P>::check_structure().map_err(PolynomialError::InvalidBaseField)?;

        match polynomial.degree() {
            None | Some(0) => return Ok(IrreducibilityStatus::Constant),
            Some(1) => return Ok(IrreducibilityStatus::Linear),
            Some(_) => {}
        }

        let original_leading = *polynomial
            .leading_coefficient()
            .expect("non-constant polynomial has a leading coefficient");
        let monic_polynomial = polynomial.make_monic()?;
        let degree = monic_polynomial
            .degree()
            .expect("non-constant polynomial has a degree");

        for divisor_degree in 1..=(degree / 2) {
            if let Some((divisor, quotient)) =
                find_monic_divisor::<P>(&monic_polynomial, &original_leading, divisor_degree)?
            {
                return Ok(IrreducibilityStatus::Reducible { divisor, quotient });
            }
        }

        Ok(IrreducibilityStatus::Irreducible)
    }
}

fn find_monic_divisor<const P: u64>(
    monic_polynomial: &DensePolynomial<Fp<P>>,
    original_leading: &<Fp<P> as Field>::Elem,
    divisor_degree: usize,
) -> Result<Option<DenseIrreducibilityWitness<P>>, PolynomialError> {
    let mut coefficients = Vec::with_capacity(divisor_degree);
    search_monic_divisor_coefficients::<P>(
        monic_polynomial,
        original_leading,
        divisor_degree,
        0,
        &mut coefficients,
    )
}

fn search_monic_divisor_coefficients<const P: u64>(
    monic_polynomial: &DensePolynomial<Fp<P>>,
    original_leading: &<Fp<P> as Field>::Elem,
    divisor_degree: usize,
    next_degree: usize,
    coefficients: &mut Vec<<Fp<P> as Field>::Elem>,
) -> Result<Option<DenseIrreducibilityWitness<P>>, PolynomialError> {
    if next_degree == divisor_degree {
        let mut divisor_coefficients = coefficients.clone();
        divisor_coefficients.push(Fp::<P>::one());
        let divisor = DensePolynomial::<Fp<P>>::new(divisor_coefficients);
        let (quotient_for_monic_polynomial, remainder) = monic_polynomial.div_rem(&divisor)?;

        if remainder.is_zero() {
            let quotient = quotient_for_monic_polynomial.scale(original_leading);
            return Ok(Some((divisor, quotient)));
        }

        return Ok(None);
    }

    for value in 0..P {
        coefficients.push(Fp::<P>::elem_from_u64(value));

        if let Some(witness) = search_monic_divisor_coefficients::<P>(
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
