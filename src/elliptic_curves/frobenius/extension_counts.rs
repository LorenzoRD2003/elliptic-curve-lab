use core::num::NonZeroU32;

use num_bigint::{BigInt, BigUint, Sign};

use crate::elliptic_curves::CurveError;
use crate::elliptic_curves::frobenius::FrobeniusTrace;
use crate::elliptic_curves::traits::EnumerableCurveModel;
use crate::fields::{EnumerableFiniteField, Field, FiniteField, FiniteFieldDescriptor, SqrtField};
use crate::numerics::OrderTwoLinearRecurrence;

/// Exact point-count report for `E(F_{q^n})` derived from the Frobenius trace.
///
/// Let `E` be an elliptic curve over `F_q`, and let `t` be the trace of the
/// relative Frobenius `π_q`. If `α, β` are the roots of `T^2 - tT + q`,
/// then `#E(F_{q^n}) = q^n + 1 - α^n - β^n`.
///
/// Writing `s_n = α^n + β^n`, we recover theorder-2 recurrence
/// `s_0 = 2`, `s_1 = t`, `s_n = t s_{n-1} - q s_{n-2}`.
///
/// This report stores the requested extension degree `n`, the exact value
/// `q^n`, the derived power sum `s_n`, and the cardinality `#E(F_{q^n})`.
///
/// The current implementation uses arbitrary-precision integers throughout, so
/// this layer does not introduce a fixed-width overflow boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusExtensionCountReport {
    frobenius_trace: FrobeniusTrace,
    extension_degree: NonZeroU32,
    extension_field_order: BigUint,
    power_sum: BigInt,
    curve_order: BigUint,
}

impl FrobeniusExtensionCountReport {
    /// Returns the base Frobenius-trace package over `F_q`.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the extension degree `n`.
    pub fn extension_degree(&self) -> NonZeroU32 {
        self.extension_degree
    }

    /// Returns the extension-field cardinality `q^n`.
    pub fn extension_field_order(&self) -> &BigUint {
        &self.extension_field_order
    }

    /// Returns the power sum `s_n = α^n + β^n`.
    pub fn power_sum(&self) -> &BigInt {
        &self.power_sum
    }

    /// Returns the derived point count `#E(F_{q^n})`.
    pub fn curve_order(&self) -> &BigUint {
        &self.curve_order
    }
}

/// Batched point-count report for
///
/// `#E(F_q), #E(F_{q^2}), ..., #E(F_{q^N})`.
///
/// This prefix report is the natural companion to
/// [`FrobeniusExtensionCountReport`]: instead of asking for one isolated
/// extension degree, it stores the whole sequence through one maximum bound.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusExtensionCountSequenceReport {
    frobenius_trace: FrobeniusTrace,
    reports: Vec<FrobeniusExtensionCountReport>,
}

impl FrobeniusExtensionCountSequenceReport {
    /// Returns the common base Frobenius-trace package over `F_q`.
    pub fn frobenius_trace(&self) -> &FrobeniusTrace {
        &self.frobenius_trace
    }

    /// Returns the stored reports for degrees `1, 2, ..., N`.
    pub fn reports(&self) -> &[FrobeniusExtensionCountReport] {
        &self.reports
    }
}

/// Comparison between a Frobenius-derived extension count and an exhaustive
/// enumerated extension count on a small represented field.
///
/// This report is intentionally pedagogical: it keeps both algorithmic paths
/// visible instead of hiding the exhaustive route behind the derived one.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrobeniusExtensionEnumerationComparisonReport {
    trace_base_field: FiniteFieldDescriptor,
    curve_base_field: FiniteFieldDescriptor,
    relative_extension_degree: NonZeroU32,
    frobenius_count: FrobeniusExtensionCountReport,
    exhaustive_curve_order: BigUint,
    agrees: bool,
}

impl FrobeniusExtensionEnumerationComparisonReport {
    /// Returns the base field `F_q` attached to the Frobenius trace.
    pub fn trace_base_field(&self) -> &FiniteFieldDescriptor {
        &self.trace_base_field
    }

    /// Returns the actually represented base field of the enumerated curve.
    pub fn curve_base_field(&self) -> &FiniteFieldDescriptor {
        &self.curve_base_field
    }

    /// Returns the relative degree `n` such that the compared curve lives over `F_{q^n}`.
    pub fn relative_extension_degree(&self) -> NonZeroU32 {
        self.relative_extension_degree
    }

    /// Returns the count derived from Frobenius data.
    pub fn frobenius_count(&self) -> &FrobeniusExtensionCountReport {
        &self.frobenius_count
    }

    /// Returns the count obtained by direct exhaustive enumeration.
    pub fn exhaustive_curve_order(&self) -> &BigUint {
        &self.exhaustive_curve_order
    }

    /// Returns whether the two routes agree.
    pub fn agrees(&self) -> bool {
        self.agrees
    }
}

impl FrobeniusTrace {
    /// Computes `#E(F_{q^n})` from the stored Frobenius trace over `F_q`.
    ///
    /// Complexity: `Θ(log n)`
    pub fn curve_order_over_extension(
        &self,
        extension_degree: NonZeroU32,
    ) -> FrobeniusExtensionCountReport {
        let recurrence = self.extension_count_recurrence();
        let extension_field_order = self.extension_field_order(extension_degree);
        let power_sum = recurrence.nth_term(u64::from(extension_degree.get()));
        let curve_order = curve_order_from_power_sum(&extension_field_order, &power_sum);

        FrobeniusExtensionCountReport {
            frobenius_trace: self.clone(),
            extension_degree,
            extension_field_order,
            power_sum,
            curve_order,
        }
    }

    /// Computes the whole prefix `#E(F_q), #E(F_{q^2}), ..., #E(F_{q^N})`
    /// from the stored Frobenius trace over `F_q`.
    ///
    /// This is the preferred surface when a caller needs many consecutive
    /// extension counts, since it advances the recurrence only once.
    ///
    /// Complexity: `Θ(N)`
    pub fn curve_orders_over_extensions_through(
        &self,
        max_extension_degree: NonZeroU32,
    ) -> FrobeniusExtensionCountSequenceReport {
        let max_degree = max_extension_degree.get() as usize;
        let recurrence = self.extension_count_recurrence();
        let power_sums = recurrence.terms_through(max_degree);
        let field_order = self.base_field_order_biguint();

        let mut reports = Vec::with_capacity(max_degree);
        let mut extension_field_order = BigUint::from(1u8);

        for degree in 1..=max_extension_degree.get() {
            extension_field_order *= &field_order;
            let extension_degree = NonZeroU32::new(degree)
                .expect("the sequence iterator only visits positive degrees");
            let power_sum = power_sums[degree as usize].clone();
            let curve_order = curve_order_from_power_sum(&extension_field_order, &power_sum);

            reports.push(FrobeniusExtensionCountReport {
                frobenius_trace: self.clone(),
                extension_degree,
                extension_field_order: extension_field_order.clone(),
                power_sum,
                curve_order,
            });
        }

        FrobeniusExtensionCountSequenceReport {
            frobenius_trace: self.clone(),
            reports,
        }
    }

    fn extension_count_recurrence(&self) -> OrderTwoLinearRecurrence<BigInt> {
        let field_order = BigInt::from(self.base_field_order_biguint());
        OrderTwoLinearRecurrence::new(
            BigInt::from(2u8),
            BigInt::from(self.trace()),
            BigInt::from(self.trace()),
            -field_order,
        )
    }

    fn extension_field_order(&self, extension_degree: NonZeroU32) -> BigUint {
        self.base_field_order_biguint().pow(extension_degree.get())
    }

    fn base_field_order_biguint(&self) -> BigUint {
        BigUint::from(self.base_field().characteristic)
            .pow(self.base_field().extension_degree.get())
    }
}

fn curve_order_from_power_sum(extension_field_order: &BigUint, power_sum: &BigInt) -> BigUint {
    biguint_from_bigint(
        &(BigInt::from(extension_field_order.clone()) + BigInt::from(1u8) - power_sum),
    )
}

/// Compares a Frobenius-derived count against direct exhaustive enumeration on
/// a small represented extension field.
///
/// Suppose `frobenius_trace` comes from a curve over `F_q`, and `curve` is the
/// same geometric model represented over an enumerable finite field `F_{q^n}`.
/// This helper computes:
///
/// - the fast count `#E(F_{q^n})` derived from the Frobenius recurrence
/// - the slow count `#E(F_{q^n})` obtained by `curve.order()`
///
/// and records whether they agree.
///
/// Complexity:
/// The Frobenius-derived side costs `Θ(log n)` exact integer operations, but
/// the overall comparison is dominated by the exhaustive path `curve.order()`,
/// which is the cost of full point enumeration on the represented small field.
pub fn compare_extension_count_with_enumeration<E: EnumerableCurveModel>(
    curve: &E,
    frobenius_trace: &FrobeniusTrace,
) -> Result<FrobeniusExtensionEnumerationComparisonReport, CurveError>
where
    E::BaseField: EnumerableFiniteField<Elem = E::Elem> + SqrtField<Elem = E::Elem> + FiniteField,
    E::Point: PartialEq,
{
    let trace_base_field = frobenius_trace.base_field().clone();
    let curve_base_field = FiniteFieldDescriptor::new(
        E::BaseField::characteristic(),
        E::BaseField::extension_degree(),
    )
    .map_err(|_| CurveError::InvalidFrobeniusBaseField {
        characteristic: E::BaseField::characteristic(),
        extension_degree: E::BaseField::extension_degree().get(),
    })?;
    let relative_extension_degree =
        relative_extension_degree(&trace_base_field, &curve_base_field)?;
    let frobenius_count = frobenius_trace.curve_order_over_extension(relative_extension_degree);
    let exhaustive_curve_order = BigUint::from(curve.order() as u64);
    let agrees = frobenius_count.curve_order() == &exhaustive_curve_order;

    Ok(FrobeniusExtensionEnumerationComparisonReport {
        trace_base_field,
        curve_base_field,
        relative_extension_degree,
        frobenius_count,
        exhaustive_curve_order,
        agrees,
    })
}

fn relative_extension_degree(
    trace_base_field: &FiniteFieldDescriptor,
    curve_base_field: &FiniteFieldDescriptor,
) -> Result<NonZeroU32, CurveError> {
    if trace_base_field.characteristic != curve_base_field.characteristic {
        return Err(CurveError::IncompatibleFrobeniusTraceBaseField {
            trace_characteristic: trace_base_field.characteristic,
            trace_extension_degree: trace_base_field.extension_degree.get(),
            curve_characteristic: curve_base_field.characteristic,
            curve_extension_degree: curve_base_field.extension_degree.get(),
        });
    }

    let trace_degree = trace_base_field.extension_degree.get();
    let curve_degree = curve_base_field.extension_degree.get();
    if !curve_degree.is_multiple_of(trace_degree) {
        return Err(CurveError::IncompatibleFrobeniusTraceBaseField {
            trace_characteristic: trace_base_field.characteristic,
            trace_extension_degree: trace_degree,
            curve_characteristic: curve_base_field.characteristic,
            curve_extension_degree: curve_degree,
        });
    }

    NonZeroU32::new(curve_degree / trace_degree).ok_or(
        CurveError::IncompatibleFrobeniusTraceBaseField {
            trace_characteristic: trace_base_field.characteristic,
            trace_extension_degree: trace_degree,
            curve_characteristic: curve_base_field.characteristic,
            curve_extension_degree: curve_degree,
        },
    )
}

fn biguint_from_bigint(value: &BigInt) -> BigUint {
    let (sign, digits) = value.to_bytes_be();
    assert_ne!(sign, Sign::Minus, "expected a non-negative integer");
    BigUint::from_bytes_be(&digits)
}
