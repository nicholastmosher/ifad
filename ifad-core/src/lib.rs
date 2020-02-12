#![deny(warnings)]
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

mod ingest;
mod models;

use crate::models::{Annotation, Gene};

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum Aspect {
    #[serde(rename = "F")]
    MolecularFunction,
    #[serde(rename = "P")]
    BiologicalProcess,
    #[serde(rename = "C")]
    CellularComponent,
}

#[derive(Debug, Hash, Copy, Clone, Eq, PartialEq, Serialize, Deserialize)]
pub enum AnnotationStatus {
    KnownExperimental,
    KnownOther,
    Unknown,
    Unannotated,
}

#[derive(Debug)]
pub struct IndexElement<'a, 'b> {
    gene: &'a Gene,
    annotations: HashMap<Aspect, HashMap<AnnotationStatus, HashSet<&'b Annotation>>>,
}

#[derive(Debug)]
pub struct GeneIndex<'a, 'b> {
    map: HashMap<&'a String, IndexElement<'a, 'b>>,
}

impl GeneIndex<'_, '_> {
    pub fn from_records<'a, 'b>(genes: &'a [Gene], annotations: &'b [Annotation]) -> GeneIndex<'a, 'b> {
        let mut map = HashMap::new();

        for gene in genes {
            let annotations = HashMap::new();
            let index_element = IndexElement { gene, annotations };
            map.insert(&gene.gene_id, index_element);
        }

        for annotation in annotations {
            let gene_id = annotation.gene_names.get(0).expect("should get gene name");
            let gene_entry = map.get_mut(&*gene_id).expect("should find gene id");

            let annotations_by_aspect = gene_entry.annotations
                .entry(annotation.aspect)
                .or_insert_with(HashMap::new);

            let annotations_by_status = annotations_by_aspect
                .entry(annotation.annotation_status)
                .or_insert_with(HashSet::new);

            annotations_by_status.insert(annotation);
        }

        GeneIndex { map }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_index() {
        let genes: Vec<Gene> = vec![
            Gene {
                gene_id: "AT1G74030".to_string(),
                gene_product_type: "protein".to_string(),
            },
            Gene {
                gene_id: "AT1G74040".to_string(),
                gene_product_type: "protein".to_string(),
            },
        ];

        let annotations: Vec<Annotation> = vec![
            Annotation {
                db: "TAIR".to_string(),
                database_id: "locus:2031476".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "IEA".to_string(),
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                annotation_status: AnnotationStatus::KnownExperimental,
                gene_names: vec!["AT1G74030".to_string()],
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
            Annotation {
                db: "TAIR".to_string(),
                database_id: "locus:2031476".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "IEA".to_string(),
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                annotation_status: AnnotationStatus::KnownOther,
                gene_names: vec!["AT1G74030".to_string()],
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
            Annotation {
                db: "TAIR".to_string(),
                database_id: "locus:2031476".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "IEA".to_string(),
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                annotation_status: AnnotationStatus::Unknown,
                gene_names: vec!["AT1G74040".to_string()],
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
        ];

        let gene_index = GeneIndex::from_records(&genes, &annotations);
        dbg!(gene_index);
    }
}
