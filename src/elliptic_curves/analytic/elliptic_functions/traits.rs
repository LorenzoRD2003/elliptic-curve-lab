use num_complex::Complex64;

use crate::elliptic_curves::analytic::elliptic_functions::EllipticFunctionTruncation;

/// Common report surface for truncated evaluations of classical elliptic
/// functions attached to a lattice `Λ`.
///
/// This trait keeps the shared metadata explicit:
/// - the input point `z`
/// - the approximate complex value
/// - the truncation policy used
/// - the number of nonzero lattice terms that were actually summed
pub trait EllipticFunctionApproximation {
    /// Returns the original evaluation point supplied by the caller.
    fn z(&self) -> &Complex64;

    /// Returns the approximate complex value produced by the truncation.
    fn value(&self) -> &Complex64;

    /// Returns the truncation policy used for this approximation.
    fn truncation(&self) -> EllipticFunctionTruncation;

    /// Returns the number of nonzero lattice terms that were summed.
    fn terms_used(&self) -> usize;
}

/// Extra capability for elliptic-function reports that track how close the
/// evaluation point came to the nearest inspected lattice pole.
pub trait HasPoleDistance {
    /// Returns the smallest Euclidean distance from the reduced evaluation
    /// point to the lattice poles inspected during the truncated evaluation.
    fn pole_distance(&self) -> f64;
}

macro_rules! impl_elliptic_function_approximation {
    ($ty:ty) => {
        impl $crate::elliptic_curves::analytic::elliptic_functions::EllipticFunctionApproximation
            for $ty
        {
            fn z(&self) -> &Complex64 {
                &self.z
            }

            fn value(&self) -> &Complex64 {
                &self.value
            }

            fn truncation(&self) -> EllipticFunctionTruncation {
                self.truncation
            }

            fn terms_used(&self) -> usize {
                self.terms_used
            }
        }
    };
}

pub(crate) use impl_elliptic_function_approximation;

macro_rules! impl_has_pole_distance {
    ($ty:ty) => {
        impl $crate::elliptic_curves::analytic::elliptic_functions::HasPoleDistance for $ty {
            fn pole_distance(&self) -> f64 {
                self.pole_distance
            }
        }
    };
}

pub(crate) use impl_has_pole_distance;
