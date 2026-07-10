use crate::isogenies::graphs::endomorphisms::CraterReport;
use crate::visualization::Visualizable;

impl Visualizable for CraterReport {
    fn format_compact(&self) -> String {
        format!(
            "crater at ℓ = {} with {} certified node(s)",
            self.prime(),
            self.nodes().len()
        )
    }

    fn describe(&self) -> String {
        [
            "Crater evidence".to_string(),
            "---------------".to_string(),
            format!("crater prime: {}", self.prime()),
            format!("certified crater nodes: {}", self.nodes().len()),
            format!("certified crater shape: {:?}", self.shape()),
            format!(
                "certified horizontal cycles: {}",
                self.horizontal_cycle_count()
            ),
        ]
        .join("\n")
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use crate::elliptic_curves::ShortWeierstrassCurve;
    use crate::fields::Fp7;
    use crate::isogenies::graphs::IsogenyGraphBuilder;
    use crate::visualization::Visualizable;

    #[test]
    fn crater_report_explanation_summarizes_certified_evidence() {
        let curve = ShortWeierstrassCurve::<Fp7>::new(Fp7::from_i64(2), Fp7::from_i64(3))
            .expect("valid F_7 curve");
        let graph = IsogenyGraphBuilder::new(curve, 3)
            .max_depth(2)
            .deduplicate_by_base_field_isomorphism(true)
            .build()
            .expect("small F_7 degree-three graph should build");
        let crater = graph
            .volcano_crater_report(&BigUint::from(3u8))
            .expect("crater report should build");

        let explanation = crater.describe();

        assert!(explanation.contains("Crater evidence"));
        assert!(explanation.contains("crater prime: 3"));
        assert!(explanation.contains("certified crater nodes: 2"));
        assert!(explanation.contains("certified crater shape: TwoVertex"));
        assert!(explanation.contains("certified horizontal cycles: 1"));
        assert!(crater.format_compact().contains("crater at ℓ = 3"));
    }
}
