use crate::{Gene, Annotation, Aspect, AnnotationStatus};
use crate::index::{GeneKey, AnnoKey, Index};
use std::collections::HashSet;
use std::convert::TryFrom;
use itertools::Itertools;
use std::borrow::Borrow;

pub struct QueryResult<IndexRef>
    where IndexRef: Borrow<Index> + Clone,
{
    index: IndexRef,
    queried_genes: HashSet<GeneKey>,
    queried_annos: HashSet<AnnoKey>,
}

impl<IndexRef> QueryResult<IndexRef>
    where IndexRef: Borrow<Index> + Clone,
{
    pub fn empty(index: IndexRef) -> QueryResult<IndexRef> {
        QueryResult {
            index,
            queried_genes: HashSet::new(),
            queried_annos: HashSet::new(),
        }
    }

    pub fn iter_genes(&self) -> impl Iterator<Item=Gene> {
        QueryResultGeneIter {
            index: self.index.clone(),
            iter: self.queried_genes.clone().into_iter(),
        }
    }

    pub fn iter_annotations(&self) -> impl Iterator<Item=Annotation> {
        QueryResultAnnotationIter {
            index: self.index.clone(),
            iter: self.queried_annos.clone().into_iter(),
        }
    }
}

struct QueryResultGeneIter<IndexRef, GKI>
    where IndexRef: Borrow<Index>,
          GKI: Iterator<Item=GeneKey>,
{
    index: IndexRef,
    iter: GKI,
}

impl<IndexRef, GKI> Iterator for QueryResultGeneIter<IndexRef, GKI>
    where IndexRef: Borrow<Index>,
          GKI: Iterator<Item=GeneKey>,
{
    type Item = Gene;

    fn next(&mut self) -> Option<Self::Item> {
        let gene_key = self.iter.next()?;
        let gene = self.index.borrow().get_gene(&gene_key);
        Some(gene.unwrap().clone())
    }
}

struct QueryResultAnnotationIter<IndexRef, AKI>
    where IndexRef: Borrow<Index>,
          AKI: Iterator<Item=AnnoKey>,
{
    index: IndexRef,
    iter: AKI,
}

impl<IndexRef, AKI> Iterator for QueryResultAnnotationIter<IndexRef, AKI>
    where IndexRef: Borrow<Index>,
          AKI: Iterator<Item=AnnoKey>,
{
    type Item = Annotation;

    fn next(&mut self) -> Option<Self::Item> {
        let anno_key = self.iter.next()?;
        let anno = self.index.borrow().get_annotation(&anno_key);
        Some(anno.unwrap().clone())
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Segment {
    aspect: Aspect,
    annotation_status: AnnotationStatus,
}

impl TryFrom<(&str, &str)> for Segment {
    type Error = ();

    fn try_from((aspect, status): (&str, &str)) -> Result<Self, Self::Error> {
        let aspect = Aspect::try_from(aspect)?;
        let status = AnnotationStatus::try_from(status)?;
        Ok(Segment { aspect, annotation_status: status })
    }
}

impl Segment {
    pub fn new(aspect: Aspect, annotation_status: AnnotationStatus) -> Self {
        Segment { aspect, annotation_status }
    }

    pub fn query<IndexRef>(&self, index: IndexRef) -> QueryResult<IndexRef>
        where IndexRef: Borrow<Index> + Clone,
    {
        // Find all genes belonging to this segment
        let queried_genes: HashSet<GeneKey> = index.borrow().gene_index
            .get(&self.aspect)
            .and_then(|statuses| statuses.get(&self.annotation_status))
            .map(IntoIterator::into_iter).into_iter()
            .flatten()
            .copied()
            .collect();

        let queried_annos: HashSet<AnnoKey> = queried_genes.iter()
            .filter_map(|gene_key| {
                index.borrow().get_gene(gene_key)
                    .and_then(|gene| {
                        let gene_id = &gene.gene_id();
                        index.borrow().anno_index.get(*gene_id)
                    })
            })
            .flat_map(|(_, annos)| annos.iter())
            .filter(|anno_key| {
                index.borrow().get_annotation(anno_key)
                    .map(|anno| anno.aspect == self.aspect
                        && anno.annotation_status == self.annotation_status)
                    .unwrap_or(false)
            })
            .copied()
            .collect();

        QueryResult {
            index,
            queried_genes,
            queried_annos,
        }
    }
}

#[derive(Debug)]
pub enum Query {
    All,
    Union(Vec<Segment>),
    Intersection(Vec<Segment>),
}

impl Query {
    pub fn execute<IndexRef>(&self, index: IndexRef) -> QueryResult<IndexRef>
        where IndexRef: Borrow<Index> + Clone,
    {
        match self {
            Query::All => query_all(index),
            Query::Union(segments) => segments.iter()
                .map(|segment| segment.query(index.clone()))
                .fold1(|a, b| union(index.clone(), a, b))
                .unwrap_or_else(|| QueryResult::empty(index)),
            Query::Intersection(segments) => segments.iter()
                .map(|segment| segment.query(index.clone()))
                .fold1(|a, b| intersect(index.clone(), a, b))
                .unwrap_or_else(|| QueryResult::empty(index)),
        }
    }
}

fn query_all<IndexRef>(index: IndexRef) -> QueryResult<IndexRef>
    where IndexRef: Borrow<Index> + Clone,
{
    let (queried_genes, queried_annos): (HashSet<GeneKey>, HashSet<AnnoKey>) =
        index.borrow().anno_index.iter()
            .flat_map(|(_, (gene, annos))| {
                annos.iter().map(move |anno| (gene, anno))
            })
            .unzip();

    QueryResult {
        index,
        queried_genes,
        queried_annos,
    }
}

fn union<IndexRef>(
    index: IndexRef,
    first: QueryResult<IndexRef>,
    second: QueryResult<IndexRef>
) -> QueryResult<IndexRef>
    where IndexRef: Borrow<Index> + Clone,
{
    let mut queried_genes = first.queried_genes;
    let mut queried_annos = first.queried_annos;

    queried_genes.extend(second.queried_genes);
    queried_annos.extend(second.queried_annos);

    QueryResult {
        index,
        queried_genes,
        queried_annos,
    }
}

fn intersect<IndexRef>(
    index: IndexRef,
    first: QueryResult<IndexRef>,
    second: QueryResult<IndexRef>
) -> QueryResult<IndexRef>
    where IndexRef: Borrow<Index> + Clone,
{
    // Take the intersection of the first and second's queried genes
    let queried_genes: HashSet<GeneKey> = first.queried_genes.iter()
        .filter(|gene_key| second.queried_genes.contains(gene_key))
        .copied()
        .collect();

    // Take the union of the first and second's queried annotations
    let union_annotations: HashSet<AnnoKey> = first.queried_annos.iter()
        .chain(second.queried_annos.iter())
        .copied()
        .collect();

    // Keep only annotations that belong to the intersected genes
    let queried_annos: HashSet<AnnoKey> = union_annotations.into_iter()
        .filter(|anno_key| {
            index.borrow().get_annotation(anno_key)
                .and_then(|anno| anno.gene_in(&index.borrow().anno_index))
                .map(|gene_key| queried_genes.contains(&gene_key))
                .unwrap_or(false)
        })
        .collect();

    QueryResult {
        index,
        queried_genes,
        queried_annos,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AnnotationRecord, GeneRecord};

    lazy_static! {
        static ref TEST_GENE_RECORDS: Vec<GeneRecord> = vec![
            /* 0 */ GeneRecord { gene_id: "AT5G48870".to_string(), gene_product_type: "protein_coding".to_string() },
            /* 1 */ GeneRecord { gene_id: "AT1G07060".to_string(), gene_product_type: "protein_coding".to_string() },
            /* 2 */ GeneRecord { gene_id: "AT4G34200".to_string(), gene_product_type: "protein_coding".to_string() },
            /* 3 */ GeneRecord { gene_id: "AT2G34580".to_string(), gene_product_type: "protein_coding".to_string() },
            /* 4 */ GeneRecord { gene_id: "AT4G30872".to_string(), gene_product_type: "other_rna".to_string() },
        ];
        static ref TEST_GENES: Vec<Gene> = TEST_GENE_RECORDS.iter()
            .map(|record| Gene::from_record(record.clone()))
            .collect();

        static ref TEST_ANNOTATION_RECORDS: Vec<AnnotationRecord> = vec![
            // AT5G48870
            /* 00 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0046540".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(),                                  aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 01 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0005634".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "SUPERSENSITIVE TO ABA AND DROUGHT 1".to_string(), alternative_gene_name: "AT5G48870|AtSAD1|LSM5|AtLSM5|SM-like 5|AT5G48870.1|K24G6.21|K24G6_21".to_string(),                              gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2156588".to_string() },
            /* 02 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0005681".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "UniProtKB-KW:KW-0747".to_string(),                                                 aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190408".to_string(), assigned_by: "UniProt".to_string(),    annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 03 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0009536".to_string(), reference: "PMID:28887381|TAIR:Publication:501776792".to_string(),  evidence_code: "HDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190404".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 04 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0006397".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "UniProtKB-KW:KW-0507".to_string(),                                                 aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190116".to_string(), assigned_by: "UniProt".to_string(),    annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 05 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0008380".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "UniProtKB-KW:KW-0508".to_string(),                                                 aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190116".to_string(), assigned_by: "UniProt".to_string(),    annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 06 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:1990726".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(),                                  aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 07 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0009737".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),    evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 08 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0003723".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),    evidence_code: "ISS".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 09 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0009414".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),    evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 10 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0003723".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(),                                  aspect: Aspect::MolecularFunction,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 11 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(),    db_object_symbol: "SAD1".to_string(), invert: "".to_string(),         go_term: "GO:0005688".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(),                                  aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(),                           alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },

            // AT1G07060
            /* 12 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "DSB formation".to_string(),                       alternative_gene_name: "AT1G07060|ATDFO|AT1G07060.1|F10K1.23|F10K1_23".to_string(),                                                     gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2007441".to_string() },
            /* 13 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0005515".to_string(), reference: "PMID:28855712|TAIR:Publication:501776744".to_string(),  evidence_code: "IPI".to_string(), /* KnownExperimental */ additional_evidence: "UniProtKB:O23277".to_string(),                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT1G07060".to_string(),                           alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(),                                               gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190919".to_string(), assigned_by: "UniProt".to_string(),    annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 14 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0051321".to_string(), reference: "PMID:22694475|TAIR:Publication:501750129".to_string(),  evidence_code: "IGI".to_string(), /* KnownExperimental */ additional_evidence: "AGI_LocusCode:AT5G54260".to_string(),                                              aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT1G07060".to_string(),                           alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(),                                               gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170516".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 15 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(),                                    alternative_gene_name: "AT1G07060|AT1G07060.2".to_string(),                                                                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556416".to_string() },
            /* 16 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(),                                    alternative_gene_name: "AT1G07060|AT1G07060.3".to_string(),                                                                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556415".to_string() },
            /* 17 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0042138".to_string(), reference: "PMID:22694475|TAIR:Publication:501750129".to_string(),  evidence_code: "IGI".to_string(), /* KnownExperimental */ additional_evidence: "AGI_LocusCode:AT5G54260".to_string(),                                              aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT1G07060".to_string(),                           alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(),                                               gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20120718".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 18 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(),    db_object_symbol: "DFO".to_string(), invert: "".to_string(),          go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(),                                    alternative_gene_name: "AT1G07060|AT1G07060.4".to_string(),                                                                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556414".to_string() },

            // AT4G34200
            /* 19 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0051287".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "InterPro:IPR006139|InterPro:IPR006140".to_string(),                                aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "InterPro".to_string(),   annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 20 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0005886".to_string(), reference: "PMID:28887381|TAIR:Publication:501776792".to_string(),  evidence_code: "HDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190404".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 21 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009507".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 22 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0005739".to_string(), reference: "PMID:14671022|TAIR:Publication:501711651".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20050120".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 23 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009570".to_string(), reference: "PMID:20061580|TAIR:Publication:501735990".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20111011".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 24 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009555".to_string(), reference: "PMID:24058165|TAIR:Publication:501756748".to_string(),  evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20131025".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 25 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009793".to_string(), reference: "PMID:24058165|TAIR:Publication:501756748".to_string(),  evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20131025".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 26 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0004617".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000107810|RGD:61987|TAIR:locus:505006128|UniProtKB:P9WNX3".to_string(), aspect: Aspect::MolecularFunction,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 27 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009507".to_string(), reference: "PMID:18431481|TAIR:Publication:501724486".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080926".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 28 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "NOT".to_string(),      go_term: "GO:0005829".to_string(), reference: "PMID:21166475|TAIR:Publication:501741191".to_string(),  evidence_code: "RCA".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20110601".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 29 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0055114".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "UniProtKB-KW:KW-0560".to_string(),                                                 aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "UniProt".to_string(),    annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 30 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0006564".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(),          evidence_code: "IEA".to_string(), /* KnownOther        */ additional_evidence: "InterPro:IPR006236".to_string(),                                                   aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "InterPro".to_string(),   annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 31 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0004617".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(),  evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 32 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009536".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 33 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009570".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "PANTHER:PTN000107896|TAIR:locus:2090649|TAIR:locus:2124266".to_string(),           aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20181029".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 34 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0000096".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(),  evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 35 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0005524".to_string(), reference: "PMID:17137349|TAIR:Publication:501720283".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080930".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 36 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0005739".to_string(), reference: "PMID:17137349|TAIR:Publication:501720283".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080930".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 37 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009536".to_string(), reference: "PMID:22923678|TAIR:Publication:501750648".to_string(),  evidence_code: "IDA".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20120831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 38 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0006520".to_string(), reference: "PMID:21873635".to_string(),                             evidence_code: "IBA".to_string(), /* KnownOther        */ additional_evidence: "MGI:MGI:1355330|PANTHER:PTN000107810|UniProtKB:P9WNX3".to_string(),                aspect: Aspect::BiologicalProcess,  unique_gene_name: "embryo sac development arrest 9".to_string(),     alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(),                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20181029".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 39 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(),    db_object_symbol: "EDA9".to_string(), invert: "".to_string(),         go_term: "GO:0009561".to_string(), reference: "PMID:15634699|TAIR:Publication:501714637".to_string(),  evidence_code: "IMP".to_string(), /* KnownExperimental */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(),                           alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(),    gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20070307".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },

            // AT2G34580
            /* 40 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(),    db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(),    go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(),                                    alternative_gene_name: "AT2G34580|AT2G34580.1|T31E10.8|T31E10_8".to_string(),                                                           gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2062355".to_string() },
            /* 41 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(),    db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(),    go_term: "GO:0003674".to_string(), reference: "TAIR:Communication:501683652".to_string(),              evidence_code: "ND".to_string(),  /* Unknown           */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT2G34580".to_string(),                           alternative_gene_name: "AT2G34580|T31E10.8|T31E10_8".to_string(),                                                                       gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20031004".to_string(), assigned_by: "TIGR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2062356".to_string() },
            /* 42 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(),    db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(),    go_term: "GO:0008150".to_string(), reference: "TAIR:Communication:1345790".to_string(),                evidence_code: "ND".to_string(),  /* Unknown           */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT2G34580".to_string(),                           alternative_gene_name: "AT2G34580|T31E10.8|T31E10_8".to_string(),                                                                       gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20061019".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2062356".to_string() },
            /* 43 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(),    db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(),    go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),          evidence_code: "ISM".to_string(), /* KnownOther        */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(),                                    alternative_gene_name: "AT2G34580|AT2G34580.2".to_string(),                                                                             gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:4515101221".to_string() },

            // AT4G30872
            /* 44 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(),    go_term: "GO:0005575".to_string(), reference: "TAIR:Communication:1345790".to_string(),                evidence_code: "ND".to_string(),  /* Unknown           */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G30872".to_string(),                           alternative_gene_name: "AT4G30872".to_string(),                                                                                         gene_product_type: "RNA".to_string(),     taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string() },
            /* 45 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(),    go_term: "GO:0003674".to_string(), reference: "TAIR:Communication:1345790".to_string(),                evidence_code: "ND".to_string(),  /* Unknown           */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G30872".to_string(),                           alternative_gene_name: "AT4G30872".to_string(),                                                                                         gene_product_type: "RNA".to_string(),     taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string() },
            /* 46 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(),    go_term: "GO:0008150".to_string(), reference: "TAIR:Communication:1345790".to_string(),                evidence_code: "ND".to_string(),  /* Unknown           */ additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G30872".to_string(),                           alternative_gene_name: "AT4G30872".to_string(),                                                                                         gene_product_type: "RNA".to_string(),     taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(),       annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string() },
        ];

        static ref EVIDENCE_CODES: &'static [&'static str] = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
        static ref TEST_ANNOTATIONS: Vec<Annotation> = TEST_ANNOTATION_RECORDS.iter()
            .map(|record| Annotation::from_record(record.clone(), &EVIDENCE_CODES))
            .collect();
    }

    #[test]
    fn test_query_all() {
        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let result = Query::All.execute(&index);

        // All of the genes from the input should appear in the query result
        assert!(TEST_GENES.iter().enumerate().all(|(i, _)| result.queried_genes.contains(&GeneKey(i))));

        // All of the annotations from the input should appear in the query result
        assert!(TEST_ANNOTATIONS.iter().enumerate().all(|(i, _)| result.queried_annos.contains(&AnnoKey(i))));
    }

    #[test]
    fn test_query_segment_bp_exp() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let segment = Segment { aspect: BiologicalProcess, annotation_status: KnownExperimental };
        let result = segment.query(&index);

        let expected_genes_vec = vec![
            GeneKey(0),
            GeneKey(1),
            GeneKey(2),
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &result.queried_genes);

        let expected_annotations_vec = vec![
            // AT5G48870
            AnnoKey(7),
            AnnoKey(9),

            // AT1G07060
            AnnoKey(14),
            AnnoKey(17),

            // AT4G34200
            AnnoKey(24),
            AnnoKey(25),
            AnnoKey(34),
            AnnoKey(39),
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &result.queried_annos);
    }

    #[test]
    fn test_query_segment_mf_other() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let segment = Segment { aspect: MolecularFunction, annotation_status: KnownOther };
        let result = segment.query(&index);

        let expected_genes_vec = vec![
            GeneKey(0),
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &result.queried_genes);
        let expected_annotations_vec = vec![
            // AT5G48870
            AnnoKey(8),
            AnnoKey(10),
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &result.queried_annos);
    }

    #[test]
    fn test_query_union() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());

        let segment_a = Segment { aspect: BiologicalProcess, annotation_status: KnownExperimental };
        let segment_b = Segment { aspect: MolecularFunction, annotation_status: KnownOther };
        let segment_c = Segment { aspect: CellularComponent, annotation_status: KnownOther };
        let query = Query::Union(vec![segment_a, segment_b, segment_c]);
        let results = query.execute(&index);

        let expected_genes_vec = vec![
            GeneKey(0),
            GeneKey(1),
            GeneKey(2),
            GeneKey(3),
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations_vec = vec![
            // AT5G48870
            AnnoKey(7),
            AnnoKey(8),
            AnnoKey(9),
            AnnoKey(10),

            // AT1G07060
            AnnoKey(12),
            AnnoKey(14),
            AnnoKey(15),
            AnnoKey(16),
            AnnoKey(17),
            AnnoKey(18),

            // AT4G34200
            AnnoKey(24),
            AnnoKey(25),
            AnnoKey(34),
            AnnoKey(39),

            // AT2G34580
            AnnoKey(40),
            AnnoKey(43),
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &results.queried_annos);
    }

    #[test]
    fn test_query_union_unknowns() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let segment_a = Segment { aspect: BiologicalProcess, annotation_status: Unknown };
        let segment_b = Segment { aspect: MolecularFunction, annotation_status: Unknown };
        let segment_c = Segment { aspect: CellularComponent, annotation_status: Unknown };
        let query = Query::Union(vec![segment_a, segment_b, segment_c]);
        let results = query.execute(&index);

        let expected_genes_vec = vec![
            GeneKey(3),
            GeneKey(4),
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations_vec = vec![
            // AT2G34580
            AnnoKey(41),
            AnnoKey(42),

            // AT4G30872
            AnnoKey(44),
            AnnoKey(45),
            AnnoKey(46),
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &results.queried_annos);
    }

    #[test]
    fn test_query_intersection() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let segment_a = Segment { aspect: CellularComponent, annotation_status: KnownOther };
        let segment_b = Segment { aspect: MolecularFunction, annotation_status: Unknown };
        let segment_c = Segment { aspect: BiologicalProcess, annotation_status: Unknown };
        let query = Query::Intersection(vec![segment_a, segment_b, segment_c]);
        let results = query.execute(&index);

        let expected_genes_vec = vec![
            // AT2G34580
            GeneKey(3),
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations_vec = vec![
            // AT2G34580
            AnnoKey(40),
            AnnoKey(41),
            AnnoKey(42),
            AnnoKey(43),
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &results.queried_annos);
    }

    #[test]
    fn test_query_intersection_empty() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(TEST_GENES.clone(), TEST_ANNOTATIONS.clone());
        let segment_a = Segment { aspect: CellularComponent, annotation_status: KnownOther };
        let segment_b = Segment { aspect: CellularComponent, annotation_status: Unknown };
        let query = Query::Intersection(vec![segment_a, segment_b]);
        let results = query.execute(&index);

        let expected_genes = HashSet::new();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations = HashSet::new();
        assert_eq!(&expected_annotations, &results.queried_annos);
    }
}
