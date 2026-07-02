use crate::fields::traits::{
    CharacteristicTwoArtinSchreierField, EnumerableFiniteField, Field, FiniteField,
    PthRootExtraction, SqrtField,
};

use super::{GeneralWeierstrassYFiberError, GeneralWeierstrassYFiberSolver, YFiberSolveResult};

/// Internal equation for one affine fiber of `x : E -> A^1` on a general
/// Weierstrass curve:
///
/// `y^2 + uy = v`.
///
/// In characteristic different from `2`, this can be solved by completing the
/// square:
///
/// `(y + u/2)^2 = v + (u/2)^2`.
///
/// In characteristic `2`, this splits into two honest subcases:
///
/// - if `u = 0`, solve `y^2 = v`
/// - if `u != 0`, divide by `u^2` and solve the Artin-Schreier equation
///   `z^2 + z = v / u^2`, then recover `y = uz`
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct GeneralWeierstrassYFiberEquation<F: Field> {
    u: F::Elem,
    v: F::Elem,
}

impl<F: Field> GeneralWeierstrassYFiberEquation<F> {
    pub(crate) fn new(u: F::Elem, v: F::Elem) -> Self {
        Self { u, v }
    }

    /// Returns the coefficient `u` in `y^2 + uy = v`.
    pub(crate) fn u(&self) -> &F::Elem {
        &self.u
    }

    /// Returns the right-hand side `v` in `y^2 + uy = v`.
    pub(crate) fn v(&self) -> &F::Elem {
        &self.v
    }

    /// Returns the shift `u / 2` used to complete the square in odd
    /// characteristic.
    pub(crate) fn odd_characteristic_shift(
        &self,
    ) -> Result<F::Elem, GeneralWeierstrassYFiberError> {
        self.ensure_odd_characteristic()?;
        let two_inverse = F::inverse(&F::from_i64(2))?;
        Ok(F::mul(self.u(), &two_inverse))
    }

    /// Returns the completed-square right-hand side `v + (u/2)^2`.
    pub(crate) fn odd_characteristic_completed_rhs(
        &self,
    ) -> Result<F::Elem, GeneralWeierstrassYFiberError> {
        let shift = self.odd_characteristic_shift()?;
        Ok(F::add(self.v(), &F::square(&shift)))
    }

    fn ensure_odd_characteristic(&self) -> Result<(), GeneralWeierstrassYFiberError> {
        if F::has_characteristic(2) {
            return Err(GeneralWeierstrassYFiberError::UnsupportedCharacteristic {
                characteristic: F::characteristic().to_biguint(),
            });
        }
        Ok(())
    }
}

impl<F: Field + SqrtField> GeneralWeierstrassYFiberEquation<F> {
    /// Returns the one or two solutions obtained by completing the square in
    /// characteristic different from `2`.
    ///
    /// When the completed-square right-hand side is not a square, this returns
    /// `Ok(None)`.
    pub(crate) fn solve_in_odd_characteristic(&self) -> YFiberSolveResult<F> {
        let shift = self.odd_characteristic_shift()?;
        let completed_rhs = self.odd_characteristic_completed_rhs()?;
        let (left_root, right_root) = match F::sqrt_pair(&completed_rhs) {
            Some(pair) => pair,
            None => return Ok(None),
        };

        Ok(Some((
            F::sub(&left_root, &shift),
            F::sub(&right_root, &shift),
        )))
    }
}

impl<F: GeneralWeierstrassYFiberSolver> GeneralWeierstrassYFiberEquation<F> {
    /// Solves the affine fiber equation through the backend-specific route
    /// chosen by `F`.
    pub(crate) fn solve(&self) -> YFiberSolveResult<F> {
        F::solve_y_fiber_equation(self)
    }
}

impl<F: FiniteField + EnumerableFiniteField + CharacteristicTwoArtinSchreierField>
    GeneralWeierstrassYFiberEquation<F>
where
    F::Elem: PthRootExtraction,
{
    /// Returns the normalized Artin-Schreier right-hand side `v / u^2` in
    /// characteristic `2` when `u != 0`.
    pub(crate) fn characteristic_two_normalized_rhs(
        &self,
    ) -> Result<Option<F::Elem>, GeneralWeierstrassYFiberError> {
        if !F::has_characteristic(2) {
            return Err(GeneralWeierstrassYFiberError::UnsupportedCharacteristic {
                characteristic: F::characteristic().to_biguint(),
            });
        }

        if F::is_zero(self.u()) {
            return Ok(None);
        }

        Ok(Some(F::div(self.v(), &F::square(self.u()))?))
    }
}
