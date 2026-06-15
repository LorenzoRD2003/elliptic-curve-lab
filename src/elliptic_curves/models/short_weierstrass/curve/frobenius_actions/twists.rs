use super::gcd_u32;
use crate::elliptic_curves::{CurveError, ShortWeierstrassCurve};
use crate::fields::traits::FiniteField;

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Returns the `p^k`-power Frobenius twist of this curve.
    ///
    /// If `E : y^2 = x^3 + ax + b` over characteristic `p`, this returns
    /// `E^(p^k) : y^2 = x^3 + a^(p^k)x + b^(p^k)`.
    ///
    /// The implementation reduces `k` modulo the represented extension degree,
    /// since over `F_{p^r}` one already has `π_p^r = id` on field elements.
    pub fn frobenius_twist_power(&self, power: u32) -> Result<Self, CurveError> {
        Self::new(
            Self::frobenius_power_element(self.a(), power),
            Self::frobenius_power_element(self.b(), power),
        )
    }

    pub(crate) fn absolute_frobenius_preserves_curve(&self, power: u32) -> Result<(), CurveError> {
        let twist = self.frobenius_twist_power(power)?;
        if F::eq(self.a(), twist.a()) && F::eq(self.b(), twist.b()) {
            Ok(())
        } else {
            Err(CurveError::AbsoluteFrobeniusDoesNotPreserveCurve { power })
        }
    }

    pub(crate) fn absolute_frobenius_period_bound(&self, power: u32) -> u32 {
        let extension_degree = F::extension_degree().get();
        let reduced_power = power % extension_degree;
        if reduced_power == 0 {
            return 1;
        }
        extension_degree / gcd_u32(extension_degree, reduced_power)
    }

    pub(crate) fn frobenius_power_element(element: &F::Elem, power: u32) -> F::Elem {
        let extension_degree = F::extension_degree().get();
        let reduced_power = power % extension_degree;
        let mut image = element.clone();
        for _ in 0..reduced_power {
            image = F::pow(&image, F::characteristic());
        }
        image
    }
}
