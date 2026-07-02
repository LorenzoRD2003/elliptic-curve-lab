use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::traits::Isogeny;
use num_bigint::BigUint;
use num_traits::One;

/// Isogeny metadata surface that keeps the separable / inseparable degree
/// factorization explicit.
pub trait DegreeFactorizedIsogeny<Domain: CurveModel, Codomain: CurveModel>:
    Isogeny<Domain, Codomain>
{
    fn separable_degree(&self) -> BigUint;

    fn inseparable_degree(&self) -> BigUint;

    fn total_degree(&self) -> BigUint {
        self.separable_degree() * self.inseparable_degree()
    }

    fn is_purely_inseparable(&self) -> bool {
        self.separable_degree().is_one() && self.inseparable_degree() > BigUint::one()
    }

    fn is_separable_by_degree(&self) -> bool {
        self.inseparable_degree().is_one()
    }
}
