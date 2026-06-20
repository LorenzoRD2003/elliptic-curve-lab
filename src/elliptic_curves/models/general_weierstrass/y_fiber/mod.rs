mod curve;
mod equation;
mod error;
mod solvers;

pub(crate) use equation::GeneralWeierstrassYFiberEquation;
pub(crate) use error::GeneralWeierstrassYFiberError;
pub(crate) use solvers::GeneralWeierstrassYFiberSolver;

#[cfg(test)]
pub(crate) use solvers::solve_in_characteristic_two;

use crate::fields::traits::Field;

type YFiberSolutions<F> = Option<(<F as Field>::Elem, <F as Field>::Elem)>;
type YFiberSolveResult<F> = Result<YFiberSolutions<F>, GeneralWeierstrassYFiberError>;
