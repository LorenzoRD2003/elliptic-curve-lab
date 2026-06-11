use crate::elliptic_curves::traits::FiniteGroupCurveModel;
use crate::fields::{EnumerableFiniteField, Field, SqrtField};
use crate::isogenies::scalar_multiplication::{
    ScalarCharacteristicFactorization, ScalarMultiplicationIsogeny,
};
use crate::isogenies::{IsogenyError, MixedKernelDescription};

impl ScalarCharacteristicFactorization {
    /// Returns the exponent `e` in the factorization `n = p^e m`
    /// where `p` is the base-field characteristic and `gcd(m, p) = 1`.
    pub fn p_power_exponent(&self) -> u32 {
        self.p_power_exponent
    }

    /// Returns the prime-to-characteristic factor `m` in
    /// `n = p^e m`.
    pub fn separable_part(&self) -> u64 {
        self.separable_part
    }

    /// Returns the degree currently attributed to the visible reduced part.
    ///
    /// In the current implementation this is `m^2`, where `m = separable_part()`.
    pub fn separable_degree(&self) -> usize {
        self.separable_degree
    }

    /// Returns the residual characteristic-`p` degree bucket.
    ///
    /// In the current implementation this is `p^(2e)`, where `n = p^e m`.
    ///
    /// Scope note:
    /// this is the degree of the full `p`-power contribution, which the
    /// current public kernel description does not yet refine into ordinary
    /// versus supersingular subcases.
    pub fn infinitesimal_degree(&self) -> usize {
        self.infinitesimal_degree
    }

    pub(super) fn new(
        p_power_exponent: u32,
        separable_part: u64,
        separable_degree: usize,
        infinitesimal_degree: usize,
    ) -> Self {
        Self {
            p_power_exponent,
            separable_part,
            separable_degree,
            infinitesimal_degree,
        }
    }
}

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
        let characteristic = C::BaseField::characteristic();
        factor_scalar_by_characteristic(self.scalar(), characteristic)
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
            return Ok(self.kernel_points.clone());
        }

        let visible_scalar = factorization.separable_part();
        let identity = self.curve.identity();
        let visible_points = self
            .curve
            .points()
            .into_iter()
            .map(|point| -> Result<_, IsogenyError> {
                let image = self.curve.mul_scalar(&point, visible_scalar)?;
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

fn factor_scalar_by_characteristic(
    scalar: u64,
    characteristic: u64,
) -> ScalarCharacteristicFactorization {
    if characteristic == 0 || characteristic == 1 {
        return ScalarCharacteristicFactorization::new(0, scalar, square_as_usize(scalar), 1);
    }

    let mut exponent = 0u32;
    let mut separable_part = scalar;
    while separable_part % characteristic == 0 {
        separable_part /= characteristic;
        exponent += 1;
    }

    let separable_degree = square_as_usize(separable_part);
    let infinitesimal_degree = pow_u64_as_usize(characteristic, exponent.saturating_mul(2));

    ScalarCharacteristicFactorization::new(
        exponent,
        separable_part,
        separable_degree,
        infinitesimal_degree,
    )
}

fn square_as_usize(value: u64) -> usize {
    usize::try_from(u128::from(value) * u128::from(value))
        .expect("educational scalar degree should fit into usize")
}

fn pow_u64_as_usize(base: u64, exponent: u32) -> usize {
    usize::try_from(u128::from(base).pow(exponent))
        .expect("educational characteristic-power degree should fit into usize")
}
