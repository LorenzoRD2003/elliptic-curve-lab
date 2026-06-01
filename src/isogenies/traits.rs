use super::IsogenyError;
use crate::elliptic_curves::{CurveModel, FiniteGroupCurveModel};
use crate::fields::{EnumerableFiniteField, SqrtField};

/// Minimal shared interface for explicit elliptic-curve isogeny objects.
///
/// This trait is intentionally austere. An isogeny is exposed only through:
///
/// - its domain curve
/// - its codomain curve
/// - its degree
/// - point evaluation
/// - the explicit kernel points used in the small finite educational setting
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

    /// Returns the explicit kernel points used by the isogeny.
    ///
    /// This deliberately exposes the small finite representation
    /// instead of hiding the kernel behind a more opaque quotient object.
    fn kernel_points(&self) -> &[Domain::Point];
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
/// - every declared kernel point maps to the codomain identity
/// - the map respects the additive group law on every pair of domain points
/// - the explicit `kernel_points()` slice matches the full identity fiber
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
                return Err(IsogenyError::ImagePointNotOnCodomain);
            }
        }

        Ok(())
    }

    /// Exhaustively checks that every declared kernel point maps to `O`.
    fn verify_maps_kernel_to_identity(&self) -> Result<(), IsogenyError> {
        let codomain_identity = self.codomain().identity();

        for point in self.kernel_points() {
            if self.evaluate(point)? != codomain_identity {
                return Err(IsogenyError::KernelPointDoesNotMapToIdentity);
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
                    return Err(IsogenyError::HomomorphismViolation);
                }
            }
        }

        Ok(())
    }

    /// Exhaustively checks that the explicit kernel equals the full identity fiber.
    fn verify_kernel_exactness(&self) -> Result<(), IsogenyError> {
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

        if explicit_point_sets_match(self.kernel_points(), &actual_kernel) {
            Ok(())
        } else {
            Err(IsogenyError::KernelMismatch)
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
