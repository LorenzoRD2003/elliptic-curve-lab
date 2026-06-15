use std::collections::BTreeMap;

use crate::elliptic_curves::{
    ShortWeierstrassCurve,
    short_weierstrass::division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
};
use crate::fields::traits::Field;
use crate::polynomials::DensePolynomial;

impl<F: Field> ShortWeierstrassCurve<F> {
    /// Returns the odd-index division polynomial `ψ_n` as an element of `F[x]`.
    pub(crate) fn odd_division_polynomial(
        &self,
        n: usize,
    ) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
        let mut odd_cache = BTreeMap::new();
        let mut even_cache = BTreeMap::new();
        self.odd_division_polynomial_inner(n, &mut odd_cache, &mut even_cache)
    }

    /// Returns the polynomial factor `f_n(x) ∈ F[x]` for an even-index division
    /// polynomial `ψ_n = y * f_n(x)`.
    pub(crate) fn even_division_polynomial_factor(
        &self,
        n: usize,
    ) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
        let mut odd_cache = BTreeMap::new();
        let mut even_cache = BTreeMap::new();
        self.even_division_polynomial_inner(n, &mut odd_cache, &mut even_cache)
    }

    fn odd_division_polynomial_inner(
        &self,
        n: usize,
        odd_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
        even_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
    ) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
        if let Some(polynomial) = odd_cache.get(&n) {
            return Ok(polynomial.clone());
        }

        let polynomial = match n {
            0 => return Err(DivisionPolynomialError::ZeroIndex),
            1 | 3 => match self.base_division_polynomial(n)? {
                DivisionPolynomialForm::InX(polynomial) => polynomial,
                DivisionPolynomialForm::YTimes(_) => {
                    unreachable!("odd base case must lie in F[x]")
                }
            },
            _ if n.is_multiple_of(2) => {
                return Err(DivisionPolynomialError::EvenIndexRequiresYFactor { n });
            }
            _ => {
                let m = (n - 1) / 2;
                let rhs_squared = &self.to_cubic().square();

                if m.is_multiple_of(2) {
                    let even_m_plus_2 =
                        self.even_division_polynomial_inner(m + 2, odd_cache, even_cache)?;
                    let even_m = self.even_division_polynomial_inner(m, odd_cache, even_cache)?;
                    let odd_m_minus_1 =
                        self.odd_division_polynomial_inner(m - 1, odd_cache, even_cache)?;
                    let odd_m_plus_1 =
                        self.odd_division_polynomial_inner(m + 1, odd_cache, even_cache)?;

                    rhs_squared
                        .mul(&even_m_plus_2.mul(&even_m.square()))
                        .sub(&odd_m_minus_1.mul(&odd_m_plus_1.cube()))
                } else {
                    let odd_m_plus_2 =
                        self.odd_division_polynomial_inner(m + 2, odd_cache, even_cache)?;
                    let odd_m = self.odd_division_polynomial_inner(m, odd_cache, even_cache)?;
                    let even_m_minus_1 =
                        self.even_division_polynomial_inner(m - 1, odd_cache, even_cache)?;
                    let even_m_plus_1 =
                        self.even_division_polynomial_inner(m + 1, odd_cache, even_cache)?;

                    odd_m_plus_2
                        .mul(&odd_m.cube())
                        .sub(&rhs_squared.mul(&even_m_minus_1.mul(&even_m_plus_1.cube())))
                }
            }
        };

        odd_cache.insert(n, polynomial.clone());
        Ok(polynomial)
    }

    fn even_division_polynomial_inner(
        &self,
        n: usize,
        odd_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
        even_cache: &mut BTreeMap<usize, DensePolynomial<F>>,
    ) -> Result<DensePolynomial<F>, DivisionPolynomialError> {
        if let Some(polynomial) = even_cache.get(&n) {
            return Ok(polynomial.clone());
        }

        let polynomial = match n {
            0 => return Err(DivisionPolynomialError::ZeroIndex),
            2 | 4 => match self.base_division_polynomial(n)? {
                DivisionPolynomialForm::YTimes(polynomial) => polynomial,
                DivisionPolynomialForm::InX(_) => {
                    unreachable!("even base case must lie in yF[x]")
                }
            },
            _ if n % 2 == 1 => return Err(DivisionPolynomialError::UnsupportedIndex { n }),
            _ => {
                let m = n / 2;
                let inverse_two = F::inverse(&F::from_i64(2)).expect(
                    "validated short-Weierstrass curves have characteristic different from 2",
                );

                if m.is_multiple_of(2) {
                    let even_m = self.even_division_polynomial_inner(m, odd_cache, even_cache)?;
                    let even_m_plus_2 =
                        self.even_division_polynomial_inner(m + 2, odd_cache, even_cache)?;
                    let even_m_minus_2 =
                        self.even_division_polynomial_inner(m - 2, odd_cache, even_cache)?;
                    let odd_m_minus_1 =
                        self.odd_division_polynomial_inner(m - 1, odd_cache, even_cache)?;
                    let odd_m_plus_1 =
                        self.odd_division_polynomial_inner(m + 1, odd_cache, even_cache)?;

                    even_m
                        .mul(
                            &even_m_plus_2
                                .mul(&odd_m_minus_1.square())
                                .sub(&even_m_minus_2.mul(&odd_m_plus_1.square())),
                        )
                        .scale(&inverse_two)
                } else {
                    let odd_m = self.odd_division_polynomial_inner(m, odd_cache, even_cache)?;
                    let odd_m_plus_2 =
                        self.odd_division_polynomial_inner(m + 2, odd_cache, even_cache)?;
                    let odd_m_minus_2 =
                        self.odd_division_polynomial_inner(m - 2, odd_cache, even_cache)?;
                    let even_m_minus_1 =
                        self.even_division_polynomial_inner(m - 1, odd_cache, even_cache)?;
                    let even_m_plus_1 =
                        self.even_division_polynomial_inner(m + 1, odd_cache, even_cache)?;

                    odd_m
                        .mul(
                            &odd_m_plus_2
                                .mul(&even_m_minus_1.square())
                                .sub(&odd_m_minus_2.mul(&even_m_plus_1.square())),
                        )
                        .scale(&inverse_two)
                }
            }
        };

        even_cache.insert(n, polynomial.clone());
        Ok(polynomial)
    }
}
