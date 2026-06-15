use num_bigint::BigUint;

use crate::elliptic_curves::{AffinePoint, CurveError, ShortWeierstrassCurve};
use crate::fields::traits::FiniteField;

use super::{
    PointOrderFromMultipleReport,
    validation::{normalized_factorization, trusted_normalized_factorization},
};

impl<F: FiniteField> ShortWeierstrassCurve<F> {
    /// Recovers the exact order of one point from a known annihilating multiple.
    ///
    /// Input contract:
    ///
    /// - `point` must lie on the curve
    /// - `multiple = M` must satisfy `[M]P = O`
    /// - `factorization` must be the prime-power factorization `M = Π ℓᵢ^eᵢ`
    ///
    /// The algorithm first isolates the `ℓ`-primary component
    /// `Q_ℓ = [M / ℓ^e]P` for each prime power `ℓ^e`, then recovers the exact
    /// local exponent by repeated multiplication by `ℓ` alone through the
    /// internal additive-group helper.
    ///
    /// Complexity:
    /// If `M = Π ℓᵢ^eᵢ` has `r` distinct prime factors and `E = Σ eᵢ`, then
    /// after factorization validation the current implementation performs:
    ///
    /// - `Θ(r)` big-scalar group multiplications to isolate the primary
    ///   components `Q_ℓ = [M / ℓ^e]P`
    /// - `Θ(E)` further group multiplications by the single primes `ℓᵢ` to
    ///   recover the exact local exponents incrementally
    ///
    /// The normalization/validation pass also performs a sort of the supplied
    /// factors and exact integer checks that the listed prime powers multiply
    /// back to `M`. In the current implementation, each group multiplication
    /// uses the internal double-and-add `BigUint` scalar path.
    pub fn point_order_from_multiple(
        &self,
        point: &AffinePoint<F>,
        multiple: BigUint,
        factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        self.validate_point_order_from_multiple_inputs(point, &multiple)?;
        let normalized_factorization =
            normalized_factorization(&multiple, factorization)?.into_factors();
        self.recover_point_order_from_normalized_factorization(
            point,
            multiple,
            &normalized_factorization,
        )
    }

    /// Internal variant for callers that already trust the prime-power
    /// factorization of the supplied multiple.
    ///
    /// This route still checks basic structural coherence such as:
    ///
    /// - non-empty factorization
    /// - positive exponents
    /// - strictly increasing prime labels after normalization
    /// - product equal to the supplied multiple
    ///
    /// But it intentionally skips primality certification of the factor bases.
    #[allow(dead_code)]
    pub(crate) fn point_order_from_multiple_with_trusted_factorization(
        &self,
        point: &AffinePoint<F>,
        multiple: BigUint,
        factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        self.validate_point_order_from_multiple_inputs(point, &multiple)?;
        let normalized_factorization =
            trusted_normalized_factorization(&multiple, factorization)?.into_factors();
        self.recover_point_order_from_normalized_factorization(
            point,
            multiple,
            &normalized_factorization,
        )
    }
}
