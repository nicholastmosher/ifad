use std::collections::{HashMap, HashSet};
use crate::{Aspect, AnnotationStatus, Gene, Annotation};

pub type GeneIndex<'a> = HashMap<Aspect, HashMap<AnnotationStatus, HashSet<&'a Gene<'a>>>>;
pub type AnnoIndex<'a, 'b> = HashMap<String, (&'a Gene<'a>, HashSet<&'b Annotation<'b>>)>;

#[derive(Debug, Eq, PartialEq)]
pub struct Index<'a, 'b> {
    pub gene_index: GeneIndex<'a>,
    pub anno_index: AnnoIndex<'a, 'b>,
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
                None => continue, // TODO collect warnings
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

        Index { gene_index, anno_index }.index_unannotated()
    }

    fn index_unannotated(mut self) -> Self {

        let aspects: &[Aspect] = &[
            Aspect::CellularComponent,
            Aspect::MolecularFunction,
            Aspect::BiologicalProcess,
        ];

        let genes_by_aspect: HashMap<Aspect, HashSet<&Gene>> = self.gene_index.iter()
            .map(|(&aspect, by_status)| {
                let genes: HashSet<_> = by_status.iter()
                    .flat_map(|(_, genes)| genes.into_iter().map(|&gene| gene)).collect();

                (aspect, genes)
            })
            .collect();

        let genes_iter = self.anno_index.iter()
            .map(|(_, (gene, _))| gene);

        let mut unannotated_by_aspect: HashMap<Aspect, HashSet<&Gene>> = HashMap::new();
        for &gene in genes_iter {
            for aspect in aspects.iter() {
                let in_aspect = genes_by_aspect.get(aspect)
                    .map(|genes| genes.contains(gene))
                    .unwrap_or(false);
                if !in_aspect {
                    unannotated_by_aspect.entry(*aspect)
                        .or_insert_with(HashSet::new)
                        .insert(gene);
                }
            }
        }

        for (aspect, genes) in unannotated_by_aspect.into_iter() {
            self.gene_index.entry(aspect)
                .or_insert_with(HashMap::new)
                .entry(AnnotationStatus::Unannotated)
                .or_insert_with(HashSet::new)
                .extend(genes);
        }

        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AnnotationRecord, GeneRecord};

    #[test]
    fn test_create_indexes() {
        let gene_records: Vec<GeneRecord> = vec![
            GeneRecord {
                gene_id: "AT1G74030".to_string(),
                gene_product_type: "protein".to_string(),
            },
            GeneRecord {
                gene_id: "AT1G74040".to_string(),
                gene_product_type: "protein".to_string(),
            },
        ];
        let genes: Vec<Gene> = gene_records.iter()
            .map(|record| Gene::from_record(record))
            .collect();

        let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
        let annotation_records: Vec<AnnotationRecord> = vec![
            AnnotationRecord {
                db: "TAIR".to_string(),
                database_id: "locus:1111111".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "EXP".to_string(), // KNOWN_EXP
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                unique_gene_name: "AT1G74030".to_string(),
                alternative_gene_name: "".to_string(),
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
            AnnotationRecord {
                db: "TAIR".to_string(),
                database_id: "locus:2222222".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "OTHER".to_string(), // KNOWN_OTHER
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                unique_gene_name: "AT1G74030".to_string(),
                alternative_gene_name: "".to_string(),
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
            AnnotationRecord {
                db: "TAIR".to_string(),
                database_id: "locus:3333333".to_string(),
                db_object_symbol: "ENO1".to_string(),
                invert: "".to_string(),
                go_term: "GO:0000015".to_string(),
                reference: "TAIR:AnalysisReference:501756966".to_string(),
                evidence_code: "ND".to_string(), // UNKNOWN
                additional_evidence: "InterPro:IPR000941".to_string(),
                aspect: Aspect::CellularComponent,
                unique_gene_name: "AT1G74040".to_string(),
                alternative_gene_name: "".to_string(),
                gene_product_type: "protein".to_string(),
                taxon: "taxon:3702".to_string(),
                date: "20190907".to_string(),
                assigned_by: "InterPro".to_string(),
                annotation_extension: "".to_string(),
                gene_product_form_id: "TAIR:locus:2031476".to_string(),
            },
        ];
        let annotations: Vec<_> = annotation_records.iter()
            .map(|record| Annotation::from_record(record, &experimental_evidence[..]))
            .collect();

        let index = Index::new(&genes, &annotations);

        let mut gene_index: GeneIndex = HashMap::new();
        let cc = gene_index.entry(Aspect::CellularComponent)
            .or_insert_with(HashMap::new);
        cc.entry(AnnotationStatus::KnownExperimental)
            .or_insert_with(HashSet::new).insert(&genes[0]);
        cc.entry(AnnotationStatus::Unknown)
            .or_insert_with(HashSet::new).insert(&genes[1]);

        // Neither gene is annotated to BiologicalProcess
        gene_index.entry(Aspect::BiologicalProcess)
            .or_insert_with(HashMap::new)
            .entry(AnnotationStatus::Unannotated)
            .or_insert_with(HashSet::new)
            .extend(&[&genes[0], &genes[1]]);

        // Neither gene is annotated to MolecularFunction
        gene_index.entry(Aspect::MolecularFunction)
            .or_insert_with(HashMap::new)
            .entry(AnnotationStatus::Unannotated)
            .or_insert_with(HashSet::new)
            .extend(&[&genes[0], &genes[1]]);

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
        assert_eq!(expected_index, index);
    }
}
