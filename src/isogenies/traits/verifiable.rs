use crate::elliptic_curves::traits::FiniteGroupCurveModel;
use crate::fields::traits::{EnumerableFiniteField, SqrtField};
use crate::isogenies::{
    error::{IsogenyError, IsogenyVerificationError},
    kernel::{KernelDescription, ReducedKernelDescription},
    traits::Isogeny,
};

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
pub trait VerifiableIsogeny<Domain: FiniteGroupCurveModel, Codomain: FiniteGroupCurveModel>:
    Isogeny<Domain, Codomain>
where
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Codomain::BaseField:
        EnumerableFiniteField<Elem = Codomain::Elem> + SqrtField<Elem = Codomain::Elem>,
    Codomain::Point: Clone + PartialEq,
{
    /// Exhaustively checks that every domain point maps into the codomain.
    fn verify_maps_domain_to_codomain(&self) -> Result<(), IsogenyError> {
        for point in self.domain().points() {
            let image = self.evaluate(&point)?;
            if !self.codomain().contains(&image) {
                return Err(IsogenyVerificationError::ImagePointNotOnCodomain.into());
            }
        }

        Ok(())
    }

    /// Exhaustively checks that every declared kernel point maps to `O`.
    fn verify_maps_kernel_to_identity(&self) -> Result<(), IsogenyError> {
        let codomain_identity = self.codomain().identity();

        for point in self.kernel_points() {
            if self.evaluate(&point)? != codomain_identity {
                return Err(IsogenyVerificationError::KernelPointDoesNotMapToIdentity.into());
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
                    return Err(IsogenyVerificationError::HomomorphismViolation.into());
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
                return Err(IsogenyVerificationError::KernelDescriptionNotPointwiseExact.into());
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
            Err(IsogenyVerificationError::KernelMismatch.into())
        }
    }
}

impl<T, Domain: FiniteGroupCurveModel, Codomain: FiniteGroupCurveModel>
    VerifiableIsogeny<Domain, Codomain> for T
where
    T: Isogeny<Domain, Codomain>,
    Domain::BaseField: EnumerableFiniteField<Elem = Domain::Elem> + SqrtField<Elem = Domain::Elem>,
    Domain::Point: Clone + PartialEq,
    Codomain::BaseField:
        EnumerableFiniteField<Elem = Codomain::Elem> + SqrtField<Elem = Codomain::Elem>,
    Codomain::Point: Clone + PartialEq,
{
}

fn explicit_point_sets_match<P>(left: &[P], right: &[P]) -> bool
where
    P: PartialEq,
{
    left.len() == right.len()
        && left.iter().all(|point| right.contains(point))
        && right.iter().all(|point| left.contains(point))
}
