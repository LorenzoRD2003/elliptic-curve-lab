use crate::elliptic_curves::traits::{CurveIsomorphism, CurveModel};
use crate::isogenies::error::{IsogenyError, IsogenyMapError, IsogenyVerificationError};

pub(crate) trait CompositionBridge<Middle: CurveModel> {
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError>;

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError>;
}

impl<Middle: CurveModel + PartialEq> CompositionBridge<Middle> for () {
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError> {
        if source != target {
            return Err(IsogenyMapError::CompositionDomainCodomainMismatch.into());
        }

        Ok(())
    }

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError> {
        Ok(point)
    }
}

impl<Middle, Bridge> CompositionBridge<Middle> for Bridge
where
    Middle: CurveModel + Clone + PartialEq,
    Bridge: CurveIsomorphism<Domain = Middle, Codomain = Middle>,
{
    fn validate_bridge(&self, source: &Middle, target: &Middle) -> Result<(), IsogenyError> {
        if self.domain() != source || self.codomain() != target {
            return Err(IsogenyMapError::CompositionDomainCodomainMismatch.into());
        }

        Ok(())
    }

    fn transport(
        &self,
        point: <Middle as CurveModel>::Point,
    ) -> Result<<Middle as CurveModel>::Point, IsogenyError> {
        self.evaluate(&point)
            .map_err(|_| IsogenyVerificationError::ImagePointNotOnCodomain.into())
    }
}
