use crate::fields::traits::*;
use core::fmt;

use crate::elliptic_curves::montgomery::{MontgomeryXzPoint, NormalizedMontgomeryCurve};

/// Error returned when one Montgomery differential-addition precondition is not
/// satisfied by the supplied `x`-line data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MontgomeryDifferentialArithmeticError {
    /// The caller supplied an identity-side differential witness that is not
    /// compatible with the requested `xADD`/`xDBLADD` special case.
    IncompatibleDifference,
}

impl fmt::Display for MontgomeryDifferentialArithmeticError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompatibleDifference => write!(
                formatter,
                "the supplied x(P-Q) data is incompatible with the requested Montgomery differential-addition special case"
            ),
        }
    }
}

impl std::error::Error for MontgomeryDifferentialArithmeticError {}

impl<F: Field> NormalizedMontgomeryCurve<F> {
    /// Returns the cached ladder coefficient
    ///
    /// `A24 = (A + 2) / 4`
    ///
    /// for the normalized Montgomery model `v^2 = x^3 + A x^2 + x`.
    fn a24(&self) -> F::Elem {
        F::div(&F::add(self.a(), &F::from_i64(2)), &F::from_i64(4))
            .expect("normalized Montgomery model lives in characteristic different from 2")
    }

    /// Doubles one projective `x`-line representative.
    ///
    /// For the normalized Montgomery model `v^2 = x^3 + A x^2 + x`, the
    /// standard `X:Z` formulas are
    ///
    /// `X([2]P) = AA * BB`, and `Z([2]P) = E * (BB + A24 * E)`,
    ///
    /// where `AA = (X + Z)^2`, `BB = (X - Z)^2`, `E = AA - BB`,
    /// and `A24 = (A + 2) / 4`.
    ///
    /// Complexity: `2` general multiplications and `2` squarings, plus a few
    /// additions/subtractions.
    pub fn x_dbl(&self, point: &MontgomeryXzPoint<F>) -> MontgomeryXzPoint<F> {
        let MontgomeryXzPoint::Finite { x, z } = point else {
            return MontgomeryXzPoint::Infinity;
        };

        let x_plus_z = F::add(x, z);
        let x_minus_z = F::sub(x, z);
        let a_sq = F::square(&x_plus_z);
        let b_sq = F::square(&x_minus_z);
        let e = F::sub(&a_sq, &b_sq);
        let x_two = F::mul(&a_sq, &b_sq);
        let z_two = F::mul(&e, &F::add(&b_sq, &F::mul(&self.a24(), &e)));

        MontgomeryXzPoint::from_xz_or_infinity(x_two, z_two)
    }

    fn x_add_identity_case(
        &self,
        left: &MontgomeryXzPoint<F>,
        right: &MontgomeryXzPoint<F>,
        difference: &MontgomeryXzPoint<F>,
    ) -> Result<Option<MontgomeryXzPoint<F>>, MontgomeryDifferentialArithmeticError> {
        match (left, right) {
            (MontgomeryXzPoint::Infinity, MontgomeryXzPoint::Infinity) => {
                Ok(Some(MontgomeryXzPoint::Infinity))
            }
            (MontgomeryXzPoint::Infinity, _) => {
                if difference == right {
                    Ok(Some(right.clone()))
                } else {
                    Err(MontgomeryDifferentialArithmeticError::IncompatibleDifference)
                }
            }
            (_, MontgomeryXzPoint::Infinity) => {
                if difference == left {
                    Ok(Some(left.clone()))
                } else {
                    Err(MontgomeryDifferentialArithmeticError::IncompatibleDifference)
                }
            }
            _ => Ok(None),
        }
    }

    fn x_add_finite(
        &self,
        left_x: &F::Elem,
        left_z: &F::Elem,
        right_x: &F::Elem,
        right_z: &F::Elem,
        difference_x: &F::Elem,
        difference_z: &F::Elem,
    ) -> MontgomeryXzPoint<F> {
        let da = F::mul(&F::sub(left_x, left_z), &F::add(right_x, right_z));
        let cb = F::mul(&F::add(left_x, left_z), &F::sub(right_x, right_z));
        let da_plus_cb = F::add(&da, &cb);
        let da_minus_cb = F::sub(&da, &cb);
        let x_sum = F::mul(difference_z, &F::square(&da_plus_cb));
        let z_sum = F::mul(difference_x, &F::square(&da_minus_cb));

        MontgomeryXzPoint::from_xz_or_infinity(x_sum, z_sum)
    }

    /// Adds two `x`-line representatives when the difference point is known.
    ///
    /// For finite `X:Z` representatives of `P`, `Q`, and `P-Q`, the formulas are
    ///
    /// `DA = (XP - ZP)(XQ + ZQ)`, `CB = (XP + ZP)(XQ - ZQ)`,
    ///
    /// `X(P+Q) = Z(P-Q) * (DA + CB)^2`, `Z(P+Q) = X(P-Q) * (DA - CB)^2`.
    ///
    /// The special cases `O + Q = Q` and `P + O = P` are handled explicitly.
    ///
    /// Complexity on the generic finite path: `6` general multiplications and
    /// `2` squarings, plus additions/subtractions.
    pub fn x_add(
        &self,
        left: &MontgomeryXzPoint<F>,
        right: &MontgomeryXzPoint<F>,
        difference: &MontgomeryXzPoint<F>,
    ) -> Result<MontgomeryXzPoint<F>, MontgomeryDifferentialArithmeticError> {
        if let Some(result) = self.x_add_identity_case(left, right, difference)? {
            return Ok(result);
        }

        let (left_x, left_z) = left
            .finite_coordinates()
            .expect("identity cases should have returned before the finite xADD path");
        let (right_x, right_z) = right
            .finite_coordinates()
            .expect("identity cases should have returned before the finite xADD path");
        let (difference_x, difference_z) = difference
            .finite_coordinates()
            .ok_or(MontgomeryDifferentialArithmeticError::IncompatibleDifference)?;

        Ok(self.x_add_finite(left_x, left_z, right_x, right_z, difference_x, difference_z))
    }

    /// Executes the combined differential step: `xDBL` on the first input and
    /// `xADD` on the pair with the supplied difference witness.
    pub fn x_dbl_add(
        &self,
        left: &MontgomeryXzPoint<F>,
        right: &MontgomeryXzPoint<F>,
        difference: &MontgomeryXzPoint<F>,
    ) -> Result<(MontgomeryXzPoint<F>, MontgomeryXzPoint<F>), MontgomeryDifferentialArithmeticError>
    {
        Ok((self.x_dbl(left), self.x_add(left, right, difference)?))
    }
}
