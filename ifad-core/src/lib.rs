#![deny(warnings)]
#![allow(dead_code)]

use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};

mod ingest;
mod models;

use crate::models::{Annotation, Gene};

pub use ingest::{AnnotationRecord, GeneRecord, MetadataReader};

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

type GeneIndex<'a> = HashMap<Aspect, HashMap<AnnotationStatus, HashSet<&'a Gene>>>;
type AnnoIndex<'a, 'b> = HashMap<String, (&'a Gene, HashSet<&'b Annotation>)>;

#[derive(Debug, Eq, PartialEq)]
pub struct Index<'a, 'b> {
    gene_index: GeneIndex<'a>,
    anno_index: AnnoIndex<'a, 'b>,
}

impl Index<'_, '_> {
    pub fn new<'a, 'b>(
        genes: &'a [Gene],
        annotations: &'b [Annotation]
    ) -> Index<'a, 'b> {
        let mut gene_index: GeneIndex = HashMap::new();
        let mut anno_index: AnnoIndex = HashMap::new();

        for gene in genes {
            anno_index.insert(gene.gene_id.to_string(), (&gene, HashSet::new()));
        }

        let mut known_other_index: GeneIndex = HashMap::new();
        for annotation in annotations {
            let gene_id = annotation.gene_in(&anno_index);
            let gene_id = match gene_id {
                Some(gene_id) => gene_id,
                None => {
                    println!("No gene for annotation {:?}", annotation);
                    continue;
                },
            };
            let (gene, gene_annotations) = anno_index
                .get_mut(&*gene_id).expect("should get gene");
            gene_annotations.insert(annotation);

            let index_to_insert =
                if annotation.annotation_status == AnnotationStatus::KnownOther {
                    &mut known_other_index
                } else {
                    &mut gene_index
                };

            index_to_insert
                .entry(annotation.aspect)
                .or_insert_with(HashMap::new)
                .entry(annotation.annotation_status)
                .or_insert_with(HashSet::new)
                .insert(&*gene);
        }

        let known_other_flat = known_other_index.into_iter()
            .flat_map(|(aspect, by_status)| {
                by_status.into_iter().flat_map(move |(_, genes)| {
                    genes.into_iter().map(move |gene| (aspect, gene))
                })
            });

        for (aspect, gene) in known_other_flat {
            let exp_for_aspect = gene_index.get(&aspect).and_then(|by_status| {
                by_status.get(&AnnotationStatus::KnownExperimental)
            });

            let exp_contains = exp_for_aspect
                .map(|exp| exp.contains(gene))
                .unwrap_or(false);

            if !exp_contains {
                gene_index.entry(aspect)
                    .or_insert_with(HashMap::new)
                    .entry(AnnotationStatus::KnownOther)
                    .or_insert_with(HashSet::new)
                    .insert(gene);
            }
        }

        Index { gene_index, anno_index }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_indexes() {
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
                database_id: "locus:1111111".to_string(),
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
                database_id: "locus:2222222".to_string(),
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
                database_id: "locus:3333333".to_string(),
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

        let index = Index::new(&genes, &annotations);

        let mut gene_index: GeneIndex = HashMap::new();
        gene_index.entry(Aspect::CellularComponent)
            .or_insert_with(HashMap::new)
            .entry(AnnotationStatus::KnownExperimental)
            .or_insert_with(HashSet::new)
            .insert(&genes[0]);
        gene_index.entry(Aspect::CellularComponent)
            .or_insert_with(HashMap::new)
            .entry(AnnotationStatus::Unknown)
            .or_insert_with(HashSet::new)
            .insert(&genes[1]);

        let mut anno_index: AnnoIndex = HashMap::new();

        let mut gene0_annotations = HashSet::new();
        gene0_annotations.insert(&annotations[0]);
        gene0_annotations.insert(&annotations[1]);
        anno_index.entry(genes[0].gene_id.to_string())
            .or_insert((&genes[0], gene0_annotations));

        let mut gene1_annotations = HashSet::new();
        gene1_annotations.insert(&annotations[2]);
        anno_index.entry(genes[1].gene_id.to_string())
            .or_insert((&genes[1], gene1_annotations));

        let expected_index = Index { gene_index, anno_index };

        dbg!(&index);
        assert_eq!(expected_index, index);
    }
}
