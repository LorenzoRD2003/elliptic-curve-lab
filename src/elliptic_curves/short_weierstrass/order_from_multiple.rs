use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::ShortWeierstrassCurve;
use crate::elliptic_curves::affine::AffinePoint;
use crate::elliptic_curves::order_from_multiple::{
    mul_scalar_biguint, recover_cyclic_primary_order,
};
use crate::elliptic_curves::traits::CurveModel;
use crate::fields::FiniteField;
use crate::numerics::{NormalizedPrimePowerFactorization, PrimePowerTable};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// One prime-by-prime reduction step in the order-from-multiple algorithm.
///
/// If the supplied multiple has the factorization `M = Π ℓᵢ^eᵢ`, the
/// short-Weierstrass wrapper isolates one prime-primary component at a time
/// and then records how much of that `ℓ`-power can be removed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointOrderReductionStep {
    prime: BigUint,
    exponent_in_multiple: u32,
    removed_exponent: u32,
    remaining_multiple_after_step: BigUint,
}

impl PointOrderReductionStep {
    /// Returns the prime `ℓ` considered at this step.
    pub fn prime(&self) -> &BigUint {
        &self.prime
    }

    /// Returns the exponent `e` of `ℓ^e` in the supplied multiple.
    pub fn exponent_in_multiple(&self) -> u32 {
        self.exponent_in_multiple
    }

    /// Returns how many copies of `ℓ` were removed while preserving
    /// annihilation of the point.
    pub fn removed_exponent(&self) -> u32 {
        self.removed_exponent
    }

    /// Returns the running remaining multiple after finishing this prime.
    pub fn remaining_multiple_after_step(&self) -> &BigUint {
        &self.remaining_multiple_after_step
    }
}

/// Report for recovering the exact order of a point from one known multiple.
///
/// Starting from one annihilating multiple `M` with `[M]P = O`, factor
///
/// `M = Π ℓᵢ^eᵢ`
///
/// and isolate each `ℓ`-primary component. The final reconstructed product of
/// the local exact powers is the exact order of `P`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PointOrderFromMultipleReport {
    supplied_multiple: BigUint,
    exact_order: BigUint,
    remaining_multiple: BigUint,
    steps: Vec<PointOrderReductionStep>,
}

impl PointOrderFromMultipleReport {
    /// Returns the original supplied multiple `M`.
    pub fn supplied_multiple(&self) -> &BigUint {
        &self.supplied_multiple
    }

    /// Returns the recovered exact order of the point.
    pub fn exact_order(&self) -> &BigUint {
        &self.exact_order
    }

    /// Returns the remaining multiple after all prime reductions.
    ///
    /// In the current algorithm this equals [`Self::exact_order`], but it is
    /// stored explicitly so the stepwise reduction story stays visible.
    pub fn remaining_multiple(&self) -> &BigUint {
        &self.remaining_multiple
    }

    /// Returns the per-prime reduction steps.
    pub fn steps(&self) -> &[PointOrderReductionStep] {
        &self.steps
    }
}

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
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        if multiple.is_zero() {
            return Err(CurveError::InvalidPointOrderMultiple {
                multiple: multiple.clone(),
            });
        }

        let normalized_factorization =
            normalized_factorization(&multiple, factorization)?.into_factors();
        ensure_annihilates(self, point, &multiple)?;

        self.point_order_from_multiple_with_normalized_factorization(
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
        if !self.contains(point) {
            return Err(CurveError::PointNotOnCurve);
        }
        if multiple.is_zero() {
            return Err(CurveError::InvalidPointOrderMultiple {
                multiple: multiple.clone(),
            });
        }

        let normalized_factorization =
            trusted_normalized_factorization(&multiple, factorization)?.into_factors();
        ensure_annihilates(self, point, &multiple)?;
        self.point_order_from_multiple_with_normalized_factorization(
            point,
            multiple,
            &normalized_factorization,
        )
    }

    fn point_order_from_multiple_with_normalized_factorization(
        &self,
        point: &AffinePoint<F>,
        supplied_multiple: BigUint,
        normalized_factorization: &[(BigUint, u32)],
    ) -> Result<PointOrderFromMultipleReport, CurveError> {
        let mut remaining_multiple = supplied_multiple.clone();
        let mut exact_order = BigUint::one();
        let mut steps = Vec::with_capacity(normalized_factorization.len());

        for (prime, exponent_in_multiple) in normalized_factorization {
            let powers = PrimePowerTable::up_through(prime, *exponent_in_multiple);
            let prime_power = powers.power(*exponent_in_multiple);
            let cofactor = &remaining_multiple / prime_power;
            let primary_component = if cofactor == BigUint::one() {
                point.clone()
            } else {
                mul_scalar_biguint(self, point, &cofactor)?
            };

            let local_report = recover_cyclic_primary_order(self, &primary_component, &powers)?;
            let removed_exponent = local_report.removed_exponent();
            let local_exact_power = powers.power(local_report.exact_exponent());
            let removed_power = powers.power(removed_exponent);

            exact_order *= local_exact_power;
            remaining_multiple /= removed_power;

            steps.push(PointOrderReductionStep {
                prime: local_report.prime().clone(),
                exponent_in_multiple: local_report.exponent_bound(),
                removed_exponent,
                remaining_multiple_after_step: remaining_multiple.clone(),
            });
        }

        Ok(PointOrderFromMultipleReport {
            supplied_multiple,
            exact_order,
            remaining_multiple,
            steps,
        })
    }
}

fn normalized_factorization(
    multiple: &BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<NormalizedPrimePowerFactorization, CurveError> {
    NormalizedPrimePowerFactorization::checked(multiple, factorization).map_err(|_| {
        CurveError::InvalidPointOrderMultipleFactorization {
            multiple: multiple.clone(),
        }
    })
}

fn trusted_normalized_factorization(
    multiple: &BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<NormalizedPrimePowerFactorization, CurveError> {
    NormalizedPrimePowerFactorization::trusted(multiple, factorization).map_err(|_| {
        CurveError::InvalidPointOrderMultipleFactorization {
            multiple: multiple.clone(),
        }
    })
}

fn ensure_annihilates<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    multiple: &BigUint,
) -> Result<(), CurveError> {
    if annihilates(curve, point, multiple)? {
        Ok(())
    } else {
        Err(CurveError::PointOrderMultipleDoesNotAnnihilatePoint {
            multiple: multiple.clone(),
        })
    }
}

fn annihilates<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    multiple: &BigUint,
) -> Result<bool, CurveError> {
    let image = mul_scalar_biguint(curve, point, multiple)?;
    Ok(curve.is_identity(&image))
}

#[cfg(test)]
pub(crate) fn point_order_from_multiple_baseline<F: FiniteField>(
    curve: &ShortWeierstrassCurve<F>,
    point: &AffinePoint<F>,
    multiple: BigUint,
    factorization: &[(BigUint, u32)],
) -> Result<PointOrderFromMultipleReport, CurveError> {
    if !curve.contains(point) {
        return Err(CurveError::PointNotOnCurve);
    }
    if multiple.is_zero() {
        return Err(CurveError::InvalidPointOrderMultiple {
            multiple: multiple.clone(),
        });
    }
    let normalized_factorization =
        normalized_factorization(&multiple, factorization)?.into_factors();
    ensure_annihilates(curve, point, &multiple)?;

    let supplied_multiple = multiple;
    let mut remaining_multiple = supplied_multiple.clone();
    let mut steps = Vec::with_capacity(normalized_factorization.len());

    for (prime, exponent_in_multiple) in &normalized_factorization {
        let mut removed_exponent = 0;

        for _ in 0..*exponent_in_multiple {
            let candidate_multiple = &remaining_multiple / prime;
            if annihilates(curve, point, &candidate_multiple)? {
                remaining_multiple = candidate_multiple;
                removed_exponent += 1;
            } else {
                break;
            }
        }

        steps.push(PointOrderReductionStep {
            prime: prime.clone(),
            exponent_in_multiple: *exponent_in_multiple,
            removed_exponent,
            remaining_multiple_after_step: remaining_multiple.clone(),
        });
    }

    Ok(PointOrderFromMultipleReport {
        supplied_multiple,
        exact_order: remaining_multiple.clone(),
        remaining_multiple,
        steps,
    })
}
