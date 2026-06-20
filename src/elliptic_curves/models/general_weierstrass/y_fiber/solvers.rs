use crate::fields::{
    ComplexApprox, FieldError, Fp, Q,
    extension_field::{ExtensionField, ExtensionFieldSpec},
    traits::{
        CharacteristicTwoArtinSchreierField, EnumerableFiniteField, Field, FiniteField,
        PthRootExtraction, SqrtField,
    },
};

use super::{GeneralWeierstrassYFiberEquation, GeneralWeierstrassYFiberError, YFiberSolveResult};

pub(crate) trait GeneralWeierstrassYFiberSolver: Field + Sized {
    fn solve_y_fiber_equation(
        equation: &GeneralWeierstrassYFiberEquation<Self>,
    ) -> YFiberSolveResult<Self>;
}

pub(crate) fn solve_in_characteristic_two<F>(
    equation: &GeneralWeierstrassYFiberEquation<F>,
) -> YFiberSolveResult<F>
where
    F: FiniteField + EnumerableFiniteField + CharacteristicTwoArtinSchreierField,
    F::Elem: PthRootExtraction,
{
    if F::characteristic() != 2 {
        return Err(GeneralWeierstrassYFiberError::UnsupportedCharacteristic {
            characteristic: F::characteristic(),
        });
    }

    if F::is_zero(equation.u()) {
        let root = equation
            .v()
            .pth_root()
            .ok_or(FieldError::Unsupported(
                "finite fields of characteristic 2 should admit unique square roots through inverse Frobenius",
            ))?;

        return Ok(Some((root.clone(), root)));
    }

    let rhs = equation
        .characteristic_two_normalized_rhs()?
        .expect("u != 0 should produce an Artin-Schreier right-hand side");
    let (left_z, right_z) = match F::solve_artin_schreier_pair(&rhs)? {
        Some(pair) => pair,
        None => return Ok(None),
    };

    let left_y = F::mul(equation.u(), &left_z);
    let right_y = F::mul(equation.u(), &right_z);
    Ok(Some((left_y, right_y)))
}

impl GeneralWeierstrassYFiberSolver for Q {
    fn solve_y_fiber_equation(
        equation: &GeneralWeierstrassYFiberEquation<Self>,
    ) -> YFiberSolveResult<Self> {
        equation.solve_in_odd_characteristic()
    }
}

impl GeneralWeierstrassYFiberSolver for ComplexApprox {
    fn solve_y_fiber_equation(
        equation: &GeneralWeierstrassYFiberEquation<Self>,
    ) -> YFiberSolveResult<Self> {
        equation.solve_in_odd_characteristic()
    }
}

impl<const P: u64> GeneralWeierstrassYFiberSolver for Fp<P>
where
    Fp<P>: SqrtField + FiniteField + EnumerableFiniteField + CharacteristicTwoArtinSchreierField,
    <Fp<P> as Field>::Elem: PthRootExtraction,
{
    fn solve_y_fiber_equation(
        equation: &GeneralWeierstrassYFiberEquation<Self>,
    ) -> YFiberSolveResult<Self> {
        if Self::characteristic() == 2 {
            solve_in_characteristic_two(equation)
        } else {
            equation.solve_in_odd_characteristic()
        }
    }
}

impl<S> GeneralWeierstrassYFiberSolver for ExtensionField<S>
where
    S: ExtensionFieldSpec,
    ExtensionField<S>:
        SqrtField + FiniteField + EnumerableFiniteField + CharacteristicTwoArtinSchreierField,
    <ExtensionField<S> as Field>::Elem: PthRootExtraction,
{
    fn solve_y_fiber_equation(
        equation: &GeneralWeierstrassYFiberEquation<Self>,
    ) -> YFiberSolveResult<Self> {
        if Self::characteristic() == 2 {
            solve_in_characteristic_two(equation)
        } else {
            equation.solve_in_odd_characteristic()
        }
    }
}
