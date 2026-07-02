use crate::elliptic_curves::traits::FiniteGroupCurveModel;
use crate::fields::traits::*;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    error::IsogenyError,
    kernel::MixedKernelDescription,
    scalar_multiplication::{ScalarCharacteristicFactorization, ScalarMultiplicationIsogeny},
};
use num_traits::ToPrimitive;

impl<C: FiniteGroupCurveModel> ScalarMultiplicationIsogeny<C>
where
    C::BaseField: EnumerableFiniteField<Elem = C::Elem> + SqrtField<Elem = C::Elem>,
    C::Point: Clone + PartialEq,
{
    /// Returns the current characteristic-based factorization data for `[n]`.
    ///
    /// Let `p` be the base-field characteristic. This method writes
    /// `n = p^e m`, with `gcd(m, p) = 1` and returns the corresponding educational
    /// bookkeeping data used by [`Self::kernel_description`].
    ///
    /// The current interpretation is:
    ///
    /// - `separable_part = m`
    /// - `separable_degree = m^2`
    /// - `infinitesimal_degree = p^(2e)`
    ///
    /// Scope note:
    /// this factorization isolates the prime-to-characteristic reduced part
    /// exactly, but it does not yet refine the characteristic-`p` contribution
    /// according to ordinary versus supersingular geometry.
    ///
    /// Complexity: `Θ(log_p n)` integer divisions plus `Θ(log e)` machine-word
    /// exponentiation.
    pub fn scalar_characteristic_factorization(&self) -> ScalarCharacteristicFactorization {
        let characteristic = C::BaseField::characteristic()
            .to_positive_biguint()
            .and_then(|value| value.to_u64())
            .expect("current scalar factorization requires characteristic to fit in u64");
        ScalarCharacteristicFactorization::from_scalar_and_characteristic(
            self.scalar(),
            characteristic,
        )
    }

    /// Returns the explicit reduced points currently visible from the
    /// prime-to-characteristic factor of `[n]`.
    ///
    /// If `n = p^e m` with `gcd(m, p) = 1`, this method enumerates the points
    /// killed by `[m]` on `E(F_q)`. In particular:
    ///
    /// - if `e = 0`, this matches the current explicit rational kernel points
    /// - if `m = 1`, the visible reduced part is just the identity point
    ///
    /// This is the reduced point data currently fed into the mixed kernel
    /// description for characteristic-divisible scalars.
    ///
    /// Complexity: `Θ(#E(F_q)) = Θ(q)` scalar-multiplication evaluations on rational
    /// points of the ambient small curve.
    pub fn visible_reduced_kernel_points(&self) -> Result<Vec<C::Point>, IsogenyError> {
        let factorization = self.scalar_characteristic_factorization();
        if factorization.p_power_exponent() == 0 {
            return Ok(self.kernel_points().to_vec());
        }

        let visible_scalar = factorization.separable_part();
        let identity = self.curve().identity();
        let visible_points = self
            .curve()
            .points()
            .into_iter()
            .map(|point| -> Result<_, IsogenyError> {
                let image = self.curve().mul_scalar(&point, visible_scalar)?;
                Ok((point, image == identity))
            })
            .collect::<Result<Vec<_>, IsogenyError>>()?
            .into_iter()
            .filter_map(|(point, kills_point)| kills_point.then_some(point))
            .collect();

        Ok(visible_points)
    }

    pub(super) fn mixed_kernel_description(
        &self,
    ) -> Result<MixedKernelDescription<C>, IsogenyError> {
        let factorization = self.scalar_characteristic_factorization();
        let reduced_points = self.visible_reduced_kernel_points()?;

        Ok(MixedKernelDescription::new(
            reduced_points,
            factorization.separable_degree(),
            factorization.infinitesimal_degree(),
            Some(format!(
                "kernel contribution from [n] = [p^{}] o [{}]",
                factorization.p_power_exponent(),
                factorization.separable_part()
            )),
        ))
    }
}
