use crate::{Aspect, AnnotationStatus, Index, Gene, Annotation};
use std::collections::HashSet;
use std::ops::Deref;

#[derive(Debug)]
pub struct QueryResult<'a> {
    genes: HashSet<&'a Gene>,
    annotations: HashSet<&'a Annotation>,
}

#[derive(Debug, Copy, Clone)]
struct Segment {
    aspect: Aspect,
    annotation_status: AnnotationStatus,
}

impl Segment {
    fn query<'a>(&self, index: &'a Index) -> QueryResult<'a> {

        // Find all genes belonging to this segment
        let genes: HashSet<&Gene> = index.gene_index
            .get(&self.aspect)
            .and_then(|statuses| statuses.get(&self.annotation_status))
            .map(IntoIterator::into_iter).into_iter()
            .flatten().map(Deref::deref)
            .collect();

        // Find all annotations belonging to those genes which
        // share the aspect and annotation of this segment
        let annotations: HashSet<&Annotation> = genes.iter()
            .map(|gene| index.anno_index.get(&gene.gene_id))
            .filter_map(|maybe_gene| maybe_gene)
            .flat_map(|(_, annos)| annos.into_iter().map(Deref::deref))
            .filter(|anno| anno.aspect == self.aspect
                && anno.annotation_status == self.annotation_status)
            .collect();

        QueryResult { genes, annotations }
    }
}

enum Query {
    All,
    Union(Vec<Segment>),
    Intersection(Vec<Segment>),
}

impl Query {
    pub fn execute<'a>(&self, _index: &'a Index) -> QueryResult<'a> {
        match self {
            Query::All => query_all(_index),
            Query::Union(_segments) => unimplemented!(),
            Query::Intersection(_segments) => unimplemented!(),
        }
    }
}

fn query_all<'a>(index: &'a Index) -> QueryResult<'a> {

    let (genes, annos): (HashSet<&Gene>, Vec<&HashSet<&Annotation>>) = index.anno_index.iter()
        .map(|(_, (gene, annos))| (gene, annos)).unzip();

    let annotations: HashSet<&Annotation> = annos.into_iter()
        .flat_map(|set| set.into_iter())
        .map(Deref::deref)
        .collect();

    QueryResult { genes, annotations }
}

fn query_union<'a>(index: &'a Index, segments: &[Segment]) -> QueryResult<'a> {

    let mut union_genes = HashSet::new();
    let mut union_annos = HashSet::new();

    for segment in segments {
        let QueryResult { genes, annotations } = segment.query(index);
        union_genes.extend(genes);
        union_annos.extend(annotations);
    }

    QueryResult { genes: union_genes, annotations: union_annos }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_all() {
        let genes = vec![
            Gene { gene_id: "AT1G74030".to_string(), gene_product_type: "protein".to_string() },
            Gene { gene_id: "AT1G74040".to_string(), gene_product_type: "protein".to_string() },
        ];

        let annotations = vec![
            Annotation { db: "TAIR".to_string(), database_id: "locus:1111111".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::KnownExperimental, gene_names: vec!["AT1G74030".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
            Annotation { db: "TAIR".to_string(), database_id: "locus:2222222".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::KnownOther, gene_names: vec!["AT1G74030".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
            Annotation { db: "TAIR".to_string(), database_id: "locus:3333333".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::Unknown, gene_names: vec!["AT1G74040".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
        ];

        let index = Index::new(&genes, &annotations);
        let result = Query::All.execute(&index);

        assert!(genes.iter().all(|gene| result.genes.contains(gene)));
        assert!(annotations.iter().all(|anno| result.annotations.contains(anno)));
    }

    #[test]
    fn test_query_segment() {
        let genes = vec![
            Gene { gene_id: "AT1G74030".to_string(), gene_product_type: "protein".to_string() },
            Gene { gene_id: "AT1G74040".to_string(), gene_product_type: "protein".to_string() },
        ];

        let annotations = vec![
            Annotation { db: "TAIR".to_string(), database_id: "locus:1111111".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::KnownExperimental, gene_names: vec!["AT1G74030".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
            Annotation { db: "TAIR".to_string(), database_id: "locus:2222222".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::KnownOther, gene_names: vec!["AT1G74030".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
            Annotation { db: "TAIR".to_string(), database_id: "locus:3333333".to_string(), db_object_symbol: "ENO1".to_string(), invert: false, go_term: "GO:0000015".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(), evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR000941".to_string(), aspect: Aspect::CellularComponent, annotation_status: AnnotationStatus::Unknown, gene_names: vec!["AT1G74040".to_string()], gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190907".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2031476".to_string() },
        ];

        let index = Index::new(&genes, &annotations);
        let segment = Segment {
            aspect: Aspect::CellularComponent,
            annotation_status: AnnotationStatus::KnownExperimental
        };
        let result = segment.query(&index);
        assert!(result.annotations.contains(&annotations[0]));
        assert!(!result.annotations.contains(&annotations[1]));
        assert!(!result.annotations.contains(&annotations[2]));
        assert!(result.genes.contains(&genes[0]));
        assert!(!result.genes.contains(&genes[1]));
    }
}
