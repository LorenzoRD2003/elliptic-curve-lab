use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::elliptic_curves::short_weierstrass::{
    division_polynomials::{DivisionPolynomialError, DivisionPolynomialForm},
    rational_torsion::{
        error::RationalTorsionError,
        integral_model::{RationalIntegralModel, integral_rational_to_bigint},
        mazur::MAZUR_CYCLIC_ORDERS,
        reduction_mod_p::small_prime_field::{ReductionPrime, ReductionResidue},
    },
};
use crate::polynomials::{IntegerPolynomial, primitive_integer_polynomial};

/// Origin of the `x`-criterion used for an order in Mazur's list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TorsionXPolynomialSource {
    /// The order-two criterion is the cubic defining the affine branch.
    TwoTorsionCubic,
    /// The criterion is the odd division polynomial `ψ_m(x)`.
    OddDivisionPolynomial,
    /// The criterion is `ψ_m/ψ_2` for an even index `m > 2`.
    ///
    /// The stored polynomial is built from the existing factor `f_m(x)` with
    /// `ψ_m = y f_m(x)`. Since `ψ_2 = 2y`, this stores the scalar-equivalent
    /// representative `2(ψ_m/ψ_2)`, which has the same roots over `ℚ` and
    /// modulo every reduction prime used here.
    EvenDivisionPolynomialOverPsi2,
}

/// Integer `x`-criterion attached to a candidate torsion order.
///
/// The stored [`TorsionXPolynomialSource`] records which mathematical
/// criterion produced the polynomial:
///
/// - order `2` uses the defining cubic `x³ + Ax + B`;
/// - odd `m` uses the division polynomial `ψ_m(x)`;
/// - even `m > 2` uses the `x`-criterion `ψ_m/ψ_2`.
///
/// In the even case the implementation obtains this criterion through the
/// division-polynomial factor already exposed by the curve layer: if
/// `ψ_m = y f_m(x)` and `ψ_2 = 2y`, then `f_m(x) = 2(ψ_m/ψ_2)`. The extra
/// scalar `2` is harmless because primitive normalization in `ℤ[x]` and root
/// testing over `ℚ` or modulo primes `p ≥ 11` are invariant under nonzero
/// scalar multiples.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct TorsionXPolynomial {
    order: usize,
    source: TorsionXPolynomialSource,
    polynomial: IntegerPolynomial,
}

impl TorsionXPolynomial {
    /// Builds the integer `x`-criterion for a Mazur-permitted non-identity order.
    ///
    /// Complexity: `Θ(n)` on the Mazur-bounded rational-torsion route, where
    /// `n` is the number of resulting coefficients cleared into `ℤ[x]`. The
    /// permitted orders form a fixed finite list, so the reused
    /// division-polynomial construction is constant-bounded here.
    pub(super) fn from_integral_model(
        model: &RationalIntegralModel,
        order: usize,
    ) -> Result<Self, TorsionXPolynomialError> {
        if !MAZUR_CYCLIC_ORDERS.contains(&order) {
            return Err(TorsionXPolynomialError::UnsupportedOrder { order });
        }

        match order {
            2 => Self::two_torsion_cubic(model),
            _ => Self::division_polynomial_x_factor(model, order),
        }
    }

    pub(super) fn order(&self) -> usize {
        self.order
    }

    pub(super) fn source(&self) -> TorsionXPolynomialSource {
        self.source
    }

    pub(super) fn polynomial(&self) -> &IntegerPolynomial {
        &self.polynomial
    }

    /// Evaluates the integer criterion modulo `p`.
    ///
    /// Complexity: `Θ(n)` modular operations for `n` dense coefficients.
    pub(super) fn evaluate_mod_prime(
        &self,
        prime: ReductionPrime,
        x: ReductionResidue,
    ) -> ReductionResidue {
        self.polynomial.to_dense_coefficients().iter().rev().fold(
            prime.zero(),
            |accumulator, coefficient| {
                let coefficient = prime.reduce_bigint(coefficient);
                prime.add(prime.mul(accumulator, x), coefficient)
            },
        )
    }

    /// Enumerates the roots of the criterion in `𝔽_p`.
    ///
    /// Complexity: `Θ(pn)` modular operations for `n` dense coefficients.
    pub(super) fn roots_mod_prime(&self, prime: ReductionPrime) -> Vec<ReductionResidue> {
        prime
            .residues()
            .filter(|x| self.evaluate_mod_prime(prime, *x).is_zero())
            .collect()
    }

    fn two_torsion_cubic(model: &RationalIntegralModel) -> Result<Self, TorsionXPolynomialError> {
        let curve = model.curve();
        let a = integral_rational_to_bigint(curve.a())?;
        let b = integral_rational_to_bigint(curve.b())?;

        Ok(Self {
            order: 2,
            source: TorsionXPolynomialSource::TwoTorsionCubic,
            polynomial: IntegerPolynomial::new(vec![b, a, BigInt::zero(), BigInt::one()]),
        })
    }

    fn division_polynomial_x_factor(
        model: &RationalIntegralModel,
        order: usize,
    ) -> Result<Self, TorsionXPolynomialError> {
        let form = model.curve().division_polynomial(order)?;
        let (source, polynomial) = match form {
            DivisionPolynomialForm::InX(polynomial) => (
                TorsionXPolynomialSource::OddDivisionPolynomial,
                primitive_integer_polynomial(&polynomial),
            ),
            DivisionPolynomialForm::YTimes(polynomial) => (
                TorsionXPolynomialSource::EvenDivisionPolynomialOverPsi2,
                primitive_integer_polynomial(&polynomial),
            ),
        };

        Ok(Self {
            order,
            source,
            polynomial,
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum TorsionXPolynomialError {
    UnsupportedOrder { order: usize },
    IntegralModel(RationalTorsionError),
    DivisionPolynomial(DivisionPolynomialError),
}

impl From<RationalTorsionError> for TorsionXPolynomialError {
    fn from(error: RationalTorsionError) -> Self {
        Self::IntegralModel(error)
    }
}

impl From<DivisionPolynomialError> for TorsionXPolynomialError {
    fn from(error: DivisionPolynomialError) -> Self {
        Self::DivisionPolynomial(error)
    }
}
