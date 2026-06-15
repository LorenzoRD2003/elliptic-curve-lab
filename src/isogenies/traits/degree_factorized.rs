use crate::elliptic_curves::traits::CurveModel;
use crate::isogenies::traits::Isogeny;

/// Isogeny metadata surface that keeps the separable / inseparable degree
/// factorization explicit.
pub trait DegreeFactorizedIsogeny<Domain: CurveModel, Codomain: CurveModel>:
    Isogeny<Domain, Codomain>
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
