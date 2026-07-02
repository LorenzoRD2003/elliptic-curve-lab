use crate::elliptic_curves::GeneralWeierstrassCurve;
use crate::fields::traits::*;

use super::{GeneralWeierstrassYFiberEquation, GeneralWeierstrassYFiberSolver, YFiberSolveResult};

impl<F: Field> GeneralWeierstrassCurve<F> {
    /// Returns the `y`-equation describing the affine fiber above one chosen
    /// `x`-coordinate:
    ///
    /// `y^2 + (a1*x + a3) y = x^3 + a2*x^2 + a4*x + a6`.
    pub(crate) fn y_fiber_equation(&self, x: &F::Elem) -> GeneralWeierstrassYFiberEquation<F> {
        let linear_coefficient = F::add(&F::mul(self.a1(), x), self.a3());
        let right_hand_side = F::add(
            &F::add(&F::cube(x), &F::mul(self.a2(), &F::square(x))),
            &F::add(&F::mul(self.a4(), x), self.a6()),
        );
        GeneralWeierstrassYFiberEquation::new(linear_coefficient, right_hand_side)
    }

    /// Solves the affine `y`-fiber above one chosen `x`-coordinate.
    pub(crate) fn solve_y_fiber(&self, x: &F::Elem) -> YFiberSolveResult<F>
    where
        F: GeneralWeierstrassYFiberSolver,
    {
        self.y_fiber_equation(x).solve()
    }
}
