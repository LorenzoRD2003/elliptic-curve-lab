use crate::elliptic_curves::traits::{CurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    IsogenyError, IsogenyVerificationError, KernelDescription, ReducedKernelDescription,
};

/// Minimal shared interface for explicit elliptic-curve isogeny objects.
///
/// This trait is intentionally austere. An isogeny is exposed only through:
///
/// - its domain curve
/// - its codomain curve
/// - its degree
/// - point evaluation
/// - an honest public description of the kernel
pub trait Isogeny<Domain, Codomain>
where
    Domain: CurveModel,
    Codomain: CurveModel,
{
    /// Returns the domain curve.
    fn domain(&self) -> &Domain;

    /// Returns the codomain curve.
    fn codomain(&self) -> &Codomain;

    /// Returns the degree of the isogeny.
    fn degree(&self) -> usize;

    /// Evaluates the isogeny at a point of the domain and returns a point of
    /// the codomain.
    ///
    /// In concrete implementations this will send every kernel point to
    /// the identity of the codomain and identify points that differ by a
    /// kernel element.
    fn evaluate(&self, point: &Domain::Point) -> Result<Codomain::Point, IsogenyError>;

    /// Returns the current public kernel description for this isogeny.
    ///
    /// This is broader than a slice of rational points: inseparable maps may
    /// have nonreduced geometric kernels that are not honestly visible through
    /// explicit point data alone.
    fn kernel_description(&self) -> KernelDescription<Domain>;

    /// Returns the explicit rational kernel points currently visible in the
    /// kernel description.
    ///
    /// This is a convenience helper for small reduced examples. It should not
    /// be read as “the full geometric kernel” in inseparable settings.
    fn kernel_points(&self) -> Vec<Domain::Point> {
        match self.kernel_description() {
            KernelDescription::Reduced(ReducedKernelDescription::RationalPointSubgroup(kernel)) => {
                kernel.points().to_vec()
            }
            KernelDescription::Reduced(
                ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints { points, .. },
            ) => points,
            KernelDescription::Mixed(description) => description.reduced_points().to_vec(),
            KernelDescription::NonReduced(_) | KernelDescription::Unknown => Vec::new(),
        }
    }
}

/// Isogeny metadata surface that keeps the separable / inseparable degree
/// factorization explicit.
pub trait DegreeFactorizedIsogeny<Domain, Codomain>: Isogeny<Domain, Codomain>
where
    Domain: CurveModel,
    Codomain: CurveModel,
{
    fn separable_degree(&self) -> u128;

    fn inseparable_degree(&self) -> u128;

    fn total_degree(&self) -> u128 {
        self.separable_degree() * self.inseparable_degree()
    }

    fn is_purely_inseparable(&self) -> bool {
        self.separable_degree() == 1 && self.inseparable_degree() > 1
    }

    fn is_separable_by_degree(&self) -> bool {
        self.inseparable_degree() == 1
    }
}

/// Exhaustive confidence checks for explicit isogenies on tiny finite curves.
///
/// This extension trait is intentionally brute-force. It is meant for the same
/// small finite educational settings as [`FiniteGroupCurveModel`], where both
/// the domain and codomain point sets can be enumerated honestly.
///
/// The default implementations verify:
///
/// - every enumerated domain point lands on the declared codomain
/// - every explicitly visible reduced kernel point maps to the codomain identity
/// - the map respects the additive group law on every pair of domain points
/// - when the kernel description is fully reduced and pointwise visible, the
///   declared kernel points match the full identity fiber
///   `{ P in E(F_q) : phi(P) = O }`
pub trait VerifiableIsogeny<Domain, Codomain>: Isogeny<Domain, Codomain>
where
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Codomain: FiniteGroupCurveModel,
    Codomain::BaseField:
        EnumerableFiniteField<Elem = Codomain::Elem> + SqrtField<Elem = Codomain::Elem>,
    Codomain::Point: Clone + PartialEq,
{
    /// Exhaustively checks that every domain point maps into the codomain.
    fn verify_maps_domain_to_codomain(&self) -> Result<(), IsogenyError> {
        for point in self.domain().points() {
            let image = self.evaluate(&point)?;
            if !self.codomain().contains(&image) {
                return Err(IsogenyError::Verification(
                    IsogenyVerificationError::ImagePointNotOnCodomain,
                ));
            }
        }

        Ok(())
    }

    /// Exhaustively checks that every declared kernel point maps to `O`.
    fn verify_maps_kernel_to_identity(&self) -> Result<(), IsogenyError> {
        let codomain_identity = self.codomain().identity();

        for point in self.kernel_points() {
            if self.evaluate(&point)? != codomain_identity {
                return Err(IsogenyError::Verification(
                    IsogenyVerificationError::KernelPointDoesNotMapToIdentity,
                ));
            }
        }

        Ok(())
    }

    /// Exhaustively checks `phi(P + Q) = phi(P) + phi(Q)` on `E(F_q)`.
    fn verify_homomorphism(&self) -> Result<(), IsogenyError> {
        for left in self.domain().points() {
            for right in self.domain().points() {
                let domain_sum = self.domain().add(&left, &right)?;
                let image_of_sum = self.evaluate(&domain_sum)?;
                let left_image = self.evaluate(&left)?;
                let right_image = self.evaluate(&right)?;
                let sum_of_images = self.codomain().add(&left_image, &right_image)?;

                if image_of_sum != sum_of_images {
                    return Err(IsogenyError::Verification(
                        IsogenyVerificationError::HomomorphismViolation,
                    ));
                }
            }
        }

        Ok(())
    }

    /// Exhaustively checks that the visible reduced kernel equals the full
    /// rational identity fiber when that comparison is mathematically honest.
    fn verify_kernel_exactness(&self) -> Result<(), IsogenyError> {
        let declared_kernel = match self.kernel_description() {
            KernelDescription::Reduced(ReducedKernelDescription::RationalPointSubgroup(kernel)) => {
                kernel.points().to_vec()
            }
            KernelDescription::Reduced(
                ReducedKernelDescription::FiniteSubgroupSchemeVisibleAsPoints { points, .. },
            ) => points,
            KernelDescription::Mixed(_)
            | KernelDescription::NonReduced(_)
            | KernelDescription::Unknown => {
                return Err(IsogenyError::Verification(
                    IsogenyVerificationError::KernelDescriptionNotPointwiseExact,
                ));
            }
        };

        let codomain_identity = self.codomain().identity();
        let actual_kernel = self
            .domain()
            .points()
            .into_iter()
            .map(|point| {
                let image = self.evaluate(&point)?;
                Ok((point, image == codomain_identity))
            })
            .collect::<Result<Vec<_>, IsogenyError>>()?
            .into_iter()
            .filter_map(|(point, maps_to_identity)| maps_to_identity.then_some(point))
            .collect::<Vec<_>>();

        if explicit_point_sets_match(declared_kernel.as_slice(), &actual_kernel) {
            Ok(())
        } else {
            Err(IsogenyError::Verification(
                IsogenyVerificationError::KernelMismatch,
            ))
        }
    }
}

impl<T, Domain, Codomain> VerifiableIsogeny<Domain, Codomain> for T
where
    T: Isogeny<Domain, Codomain>,
    Domain: FiniteGroupCurveModel,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Codomain: FiniteGroupCurveModel,
    Codomain::BaseField:
        EnumerableFiniteField<Elem = Codomain::Elem> + SqrtField<Elem = Codomain::Elem>,
    Codomain::Point: Clone + PartialEq,
{
}

fn explicit_point_sets_match<Point>(left: &[Point], right: &[Point]) -> bool
where
    Point: PartialEq,
{
    left.len() == right.len()
        && left.iter().all(|point| right.contains(point))
        && right.iter().all(|point| left.contains(point))
}

#[cfg(test)]
mod tests;
