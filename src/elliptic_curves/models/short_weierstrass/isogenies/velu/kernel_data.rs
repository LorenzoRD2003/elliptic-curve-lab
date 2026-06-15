use crate::elliptic_curves::AffinePoint;
use crate::fields::traits::Field;

#[derive(Clone, Debug)]
pub(crate) struct VeluKernelData<F: Field> {
    kernel_nonzero_points: Vec<AffinePoint<F>>,
    a_correction_sum: F::Elem,
    b_correction_sum: F::Elem,
}

impl<F: Field> VeluKernelData<F> {
    pub(crate) fn new(
        kernel_nonzero_points: Vec<AffinePoint<F>>,
        a_correction_sum: F::Elem,
        b_correction_sum: F::Elem,
    ) -> Self {
        Self {
            kernel_nonzero_points,
            a_correction_sum,
            b_correction_sum,
        }
    }

    pub(crate) fn kernel_nonzero_points(&self) -> &[AffinePoint<F>] {
        &self.kernel_nonzero_points
    }

    pub(crate) fn a_correction_sum(&self) -> &F::Elem {
        &self.a_correction_sum
    }

    pub(crate) fn b_correction_sum(&self) -> &F::Elem {
        &self.b_correction_sum
    }
}
