use std::collections::{HashMap, HashSet};
use crate::{Aspect, AnnotationStatus, Gene, Annotation};

pub type GeneIndex<'a> = HashMap<Aspect, HashMap<AnnotationStatus, HashSet<&'a Gene<'a>>>>;
pub type AnnoIndex<'a, 'b> = HashMap<String, (&'a Gene<'a>, HashSet<&'b Annotation<'b>>)>;

#[derive(Debug, Eq, PartialEq)]
pub struct Index<'a, 'b> {
    pub genes: &'a [Gene<'a>],
    pub annotations: &'b [Annotation<'b>],
    pub gene_index: GeneIndex<'a>,
    pub anno_index: AnnoIndex<'a, 'b>,
}

impl Index<'_, '_> {

    /// Create a new Index from a slice of Genes and a slice of Annotations.
    ///
    /// An Index is basically a data structure that organizes references to
    /// Genes and Annotations based on their properties. For example, references
    /// to Genes are organized by the Aspect and Annotation Status that the
    /// gene belongs to (according to the given Annotations). Additionally,
    /// sets of Annotation references are stored according to the Gene which
    /// they annotate. By organizing the data in this way, we help to optimize
    /// the speed of lookups.
    pub fn new<'a, 'b>(
        genes: &'a [Gene],
        annotations: &'b [Annotation]
    ) -> Index<'a, 'b> {
        let mut gene_index: GeneIndex = HashMap::new();
        let mut anno_index: AnnoIndex = HashMap::new();

        // The annotation index should have a key for each Gene that exists.
        // There may not exist any annotations for some genes, but we still
        // make empty entries for those genes anyways.
        for gene in genes {
            anno_index.insert(gene.gene_id.to_string(), (&gene, HashSet::new()));
        }

        // First pass: Put all KnownExperimental and Unknown annotations
        // directly into the index, but put all KnownOther annotations into
        // a temporary index.
        //
        // We will come back for a second pass to determine whether each of the
        // KnownOther annotations should be placed in the permanent index.
        let mut known_other_index: GeneIndex = HashMap::new();
        for annotation in annotations {
            let gene_id = annotation.gene_in(&anno_index)
                .map(|gene| gene.gene_id.to_string());
            let gene_id = match gene_id {
                Some(gene_id) => gene_id,
                None => continue, // TODO collect warnings
            };
            let (gene, gene_annotations) = anno_index
                .get_mut(&*gene_id).expect("should get gene");
            gene_annotations.insert(annotation);

            // Insert into temporary index for KnownOther, or
            // permanent index for KnownExperimental and Unknown.
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

        // Create an iterator over all Genes in the temporary KnownOther
        // index where each Gene is paired with the Aspect it was annotated with
        let known_other_flat = known_other_index.into_iter()
            .flat_map(|(aspect, by_status)| {
                by_status.into_iter().flat_map(move |(_, genes)| {
                    genes.into_iter().map(move |gene| (aspect, gene))
                })
            });

        // Second Pass: Adding KnownOther genes to the permanent index.
        //
        // For each Gene (G) with Aspect (A) in the KnownOther index:
        //   * Look up all KnownExperimental genes in aspect A in the permanent index
        //   * Determine whether gene G appears in that KnownExperimental set
        //   * If gene G does not appear in the KnownExperimental set, add G to
        //     the KnownOther set for aspect A in the permanent index.
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

        Index { genes, annotations, gene_index, anno_index }.index_unannotated()
    }

    /// Calculates the Unannotated section for each Aspect in the index.
    ///
    /// After an Index has been constructed with all of the Annotated categories -
    /// i.e. KnownExperimental, KnownOther, and Unknown - then we can do another
    /// pass in order to calculate the genes which are _not_ annotated to each
    /// aspect.
    fn index_unannotated(mut self) -> Self {

        let aspects: &[Aspect] = &[
            Aspect::CellularComponent,
            Aspect::MolecularFunction,
            Aspect::BiologicalProcess,
        ];

        // For each Aspect, collect a set of all Genes which are annotated to it
        let genes_by_aspect: HashMap<Aspect, HashSet<&Gene>> = self.gene_index.iter()
            .map(|(&aspect, by_status)| {
                let genes: HashSet<_> = by_status.iter()
                    .flat_map(|(_, genes)| genes.iter().copied()).collect();

                (aspect, genes)
            })
            .collect();

        // Create an iterator over _all_ genes
        let genes_iter = self.anno_index.iter()
            .map(|(_, (gene, _))| gene);

        // For each Gene (G), and for each Aspect (A):
        // If gene G does not appear in the annotations for aspect A, then
        // add gene G to the "Unannotated" set for aspect A.
        for &gene in genes_iter {
            for aspect in aspects.iter() {
                let in_aspect = genes_by_aspect.get(aspect)
                    .map(|genes| genes.contains(gene))
                    .unwrap_or(false);
                if !in_aspect {
                    self.gene_index.entry(*aspect)
                        .or_insert_with(HashMap::new)
                        .entry(AnnotationStatus::Unannotated)
                        .or_insert_with(HashSet::new)
                        .insert(gene);
                }
            }
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

        let expected_index = Index {
            genes: &genes,
            annotations: &annotations,
            gene_index,
            anno_index,
        };
        assert_eq!(expected_index, index);
    }
}
