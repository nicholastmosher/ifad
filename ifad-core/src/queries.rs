use crate::{Aspect, AnnotationStatus, Index, Gene, Annotation};
use std::collections::HashSet;
use std::ops::Deref;

#[derive(Debug)]
pub struct QueryResult<'a> {
    ordered: bool,
    genes: &'a [Gene<'a>],
    annotations: &'a [Annotation<'a>],
    queried_genes: HashSet<&'a Gene<'a>>,
    queried_annotations: HashSet<&'a Annotation<'a>>,
}

enum EitherIter<A: Iterator, B: Iterator> {
    First(A),
    Second(B),
}

impl<A, B, T> Iterator for EitherIter<A, B>
    where A: Iterator<Item=T>,
          B: Iterator<Item=T>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            EitherIter::First(a) => a.next(),
            EitherIter::Second(b) => b.next(),
        }
    }
}

impl QueryResult<'_> {
    pub fn genes_iter(&self) -> impl Iterator<Item=&Gene> {
        if self.ordered {
            EitherIter::First(self.genes.into_iter()
                .filter(move |&gene| self.queried_genes.contains(gene)))
        } else {
            EitherIter::Second(self.queried_genes.iter().map(|&gene| gene))
        }
    }

    pub fn annotations_iter(&self) -> impl Iterator<Item=&Annotation> {
        if self.ordered {
            EitherIter::First(self.annotations.into_iter()
                .filter(move |&anno| self.queried_annotations.contains(anno)))
        } else {
            EitherIter::Second(self.queried_annotations.iter().map(|&anno| anno))
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Segment {
    aspect: Aspect,
    annotation_status: AnnotationStatus,
}

impl Segment {
    pub fn new(aspect: Aspect, annotation_status: AnnotationStatus) -> Self {
        Segment { aspect, annotation_status }
    }

    pub fn query<'a>(&self, index: &'a Index) -> QueryResult<'a> {

        // Find all genes belonging to this segment
        let queried_genes: HashSet<&Gene> = index.gene_index
            .get(&self.aspect)
            .and_then(|statuses| statuses.get(&self.annotation_status))
            .map(IntoIterator::into_iter).into_iter()
            .flatten().map(Deref::deref)
            .collect();

        // Find all annotations belonging to those genes which
        // share the aspect and annotation of this segment
        let queried_annotations: HashSet<&Annotation> = queried_genes.iter()
            .map(|gene| index.anno_index.get(gene.gene_id))
            .filter_map(|maybe_gene| maybe_gene)
            .flat_map(|(_, annos)| annos.into_iter().map(Deref::deref))
            .filter(|anno| anno.aspect == self.aspect
                && anno.annotation_status == self.annotation_status)
            .collect();

        QueryResult {
            ordered: false,
            genes: &index.genes,
            annotations: &index.annotations,
            queried_genes,
            queried_annotations,
        }
    }
}

pub enum Query {
    All,
    Union(Vec<Segment>),
    Intersection(Vec<Segment>),
}

impl Query {
    pub fn execute<'a>(&self, index: &'a Index) -> QueryResult<'a> {
        match self {
            Query::All => query_all(index),
            Query::Union(segments) => query_union(index, &segments),
            Query::Intersection(_segments) => unimplemented!(),
        }
    }
}

fn query_all<'a>(index: &'a Index) -> QueryResult<'a> {

    let (queried_genes, annos): (HashSet<&Gene>, Vec<&HashSet<&Annotation>>) = index.anno_index.iter()
        .map(|(_, (gene, annos))| (gene, annos)).unzip();

    let queried_annotations: HashSet<&Annotation> = annos.into_iter()
        .flat_map(|set| set.into_iter())
        .map(Deref::deref)
        .collect();

    QueryResult {
        ordered: true,
        genes: &index.genes,
        annotations: &index.annotations,
        queried_genes,
        queried_annotations,
    }
}

fn query_union<'a>(index: &'a Index, segments: &[Segment]) -> QueryResult<'a> {
    let mut union_genes = HashSet::new();
    let mut union_annos = HashSet::new();

    for segment in segments {
        let QueryResult { queried_genes, queried_annotations, .. } = segment.query(index);
        union_genes.extend(queried_genes);
        union_annos.extend(queried_annotations);
    }

    QueryResult {
        ordered: false,
        genes: &index.genes,
        annotations: &index.annotations,
        queried_genes: union_genes,
        queried_annotations: union_annos,
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
        static ref TEST_GENES: Vec<Gene<'static>> = TEST_GENE_RECORDS.iter()
            .map(|record| Gene::from_record(record))
            .collect();

        static ref TEST_ANNOTATION_RECORDS: Vec<AnnotationRecord> = vec![
            // AT5G48870
            /* 00 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0046540".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 01 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0005634".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                                    aspect: Aspect::CellularComponent,  unique_gene_name: "SUPERSENSITIVE TO ABA AND DROUGHT 1".to_string(), alternative_gene_name: "AT5G48870|AtSAD1|LSM5|AtLSM5|SM-like 5|AT5G48870.1|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2156588".to_string() },
            /* 02 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0005681".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "UniProtKB-KW:KW-0747".to_string(),                aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190408".to_string(), assigned_by: "UniProt".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 03 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0009536".to_string(), reference: "PMID:28887381|TAIR:Publication:501776792".to_string(), evidence_code: "HDA".to_string(), additional_evidence: "".to_string(),                                    aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190404".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 04 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0006397".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "UniProtKB-KW:KW-0507".to_string(),                aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190116".to_string(), assigned_by: "UniProt".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 05 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0008380".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "UniProtKB-KW:KW-0508".to_string(),                aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190116".to_string(), assigned_by: "UniProt".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 06 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:1990726".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 07 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0009737".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),   evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                    aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 08 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0003723".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),   evidence_code: "ISS".to_string(), additional_evidence: "".to_string(),                                    aspect: Aspect::MolecularFunction,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 09 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0009414".to_string(), reference: "PMID:11740939|TAIR:Publication:1546233".to_string(),   evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                    aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20030722".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2156589".to_string() },
            /* 10 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0003723".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(), aspect: Aspect::MolecularFunction,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 11 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2156589".to_string(), db_object_symbol: "SAD1".to_string(), invert: "".to_string(), go_term: "GO:0005688".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000469347|SGD:S000000948".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "AT5G48870".to_string(), alternative_gene_name: "AT5G48870|SAD1|AtSAD1|LSM5|AtLSM5|SUPERSENSITIVE TO ABA AND DROUGHT 1|SM-like 5|K24G6.21|K24G6_21".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },

            // AT1G07060
            /* 12 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                        aspect: Aspect::CellularComponent,  unique_gene_name: "DSB formation".to_string(), alternative_gene_name: "AT1G07060|ATDFO|AT1G07060.1|F10K1.23|F10K1_23".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2007441".to_string() },
            /* 13 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0005515".to_string(), reference: "PMID:28855712|TAIR:Publication:501776744".to_string(), evidence_code: "IPI".to_string(), additional_evidence: "UniProtKB:O23277".to_string(),        aspect: Aspect::MolecularFunction,  unique_gene_name: "AT1G07060".to_string(), alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190919".to_string(), assigned_by: "UniProt".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 14 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0051321".to_string(), reference: "PMID:22694475|TAIR:Publication:501750129".to_string(), evidence_code: "IGI".to_string(), additional_evidence: "AGI_LocusCode:AT5G54260".to_string(), aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT1G07060".to_string(), alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170516".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 15 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                        aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(), alternative_gene_name: "AT1G07060|AT1G07060.2".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556416".to_string() },
            /* 16 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                        aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(), alternative_gene_name: "AT1G07060|AT1G07060.3".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556415".to_string() },
            /* 17 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0042138".to_string(), reference: "PMID:22694475|TAIR:Publication:501750129".to_string(), evidence_code: "IGI".to_string(), additional_evidence: "AGI_LocusCode:AT5G54260".to_string(), aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT1G07060".to_string(), alternative_gene_name: "AT1G07060|DFO|ATDFO|DSB formation|F10K1.23|F10K1_23".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20120718".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2007442".to_string() },
            /* 18 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2007442".to_string(), db_object_symbol: "DFO".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                        aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(), alternative_gene_name: "AT1G07060|AT1G07060.4".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:6532556414".to_string() },

            // AT4G34200
            /* 19 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0051287".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR006139|InterPro:IPR006140".to_string(),                                aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 20 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0005886".to_string(), reference: "PMID:28887381|TAIR:Publication:501776792".to_string(), evidence_code: "HDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190404".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 21 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009507".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(),         evidence_code: "ISM".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 22 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "PMID:14671022|TAIR:Publication:501711651".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20050120".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 23 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009570".to_string(), reference: "PMID:20061580|TAIR:Publication:501735990".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20111011".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 24 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009555".to_string(), reference: "PMID:24058165|TAIR:Publication:501756748".to_string(), evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20131025".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 25 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009793".to_string(), reference: "PMID:24058165|TAIR:Publication:501756748".to_string(), evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20131025".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 26 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0004617".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000107810|RGD:61987|TAIR:locus:505006128|UniProtKB:P9WNX3".to_string(), aspect: Aspect::MolecularFunction,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20170228".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 27 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009507".to_string(), reference: "PMID:18431481|TAIR:Publication:501724486".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080926".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 28 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "NOT".to_string(), go_term: "GO:0005829".to_string(), reference: "PMID:21166475|TAIR:Publication:501741191".to_string(),  evidence_code: "RCA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20110601".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 29 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0055114".to_string(), reference: "TAIR:AnalysisReference:501756968".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "UniProtKB-KW:KW-0560".to_string(),                                                 aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "UniProt".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 30 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0006564".to_string(), reference: "TAIR:AnalysisReference:501756966".to_string(),         evidence_code: "IEA".to_string(), additional_evidence: "InterPro:IPR006236".to_string(),                                                   aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190514".to_string(), assigned_by: "InterPro".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 31 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0004617".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(), evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 32 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009536".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 33 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009570".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "PANTHER:PTN000107896|TAIR:locus:2090649|TAIR:locus:2124266".to_string(),           aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20181029".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 34 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0000096".to_string(), reference: "PMID:30787133|TAIR:Publication:501784092".to_string(), evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20190510".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 35 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0005524".to_string(), reference: "PMID:17137349|TAIR:Publication:501720283".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::MolecularFunction,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080930".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 36 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "PMID:17137349|TAIR:Publication:501720283".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20080930".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:3438887".to_string() },
            /* 37 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009536".to_string(), reference: "PMID:22923678|TAIR:Publication:501750648".to_string(), evidence_code: "IDA".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20120831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },
            /* 38 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0006520".to_string(), reference: "PMID:21873635".to_string(),                            evidence_code: "IBA".to_string(), additional_evidence: "MGI:MGI:1355330|PANTHER:PTN000107810|UniProtKB:P9WNX3".to_string(),                aspect: Aspect::BiologicalProcess,  unique_gene_name: "embryo sac development arrest 9".to_string(), alternative_gene_name: "AT4G34200|PGDH1|phosphoglycerate dehydrogenase 1|AT4G34200.1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20181029".to_string(), assigned_by: "GO_Central".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "".to_string() },
            /* 39 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2124266".to_string(), db_object_symbol: "EDA9".to_string(), invert: "".to_string(), go_term: "GO:0009561".to_string(), reference: "PMID:15634699|TAIR:Publication:501714637".to_string(), evidence_code: "IMP".to_string(), additional_evidence: "".to_string(),                                                                     aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G34200".to_string(), alternative_gene_name: "AT4G34200|EDA9|PGDH1|embryo sac development arrest 9|phosphoglycerate dehydrogenase 1|F10M10.7".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20070307".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2124266".to_string() },

            // AT2G34580
            /* 40 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(), db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(), evidence_code: "ISM".to_string(), additional_evidence: "".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(), alternative_gene_name: "AT2G34580|AT2G34580.1|T31E10.8|T31E10_8".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:2062355".to_string() },
            /* 41 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(), db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(), go_term: "GO:0003674".to_string(), reference: "TAIR:Communication:501683652".to_string(),     evidence_code: "ND".to_string(),  additional_evidence: "".to_string(), aspect: Aspect::MolecularFunction,  unique_gene_name: "AT2G34580".to_string(), alternative_gene_name: "AT2G34580|T31E10.8|T31E10_8".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20031004".to_string(), assigned_by: "TIGR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2062356".to_string() },
            /* 42 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(), db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(), go_term: "GO:0008150".to_string(), reference: "TAIR:Communication:1345790".to_string(),       evidence_code: "ND".to_string(),  additional_evidence: "".to_string(), aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT2G34580".to_string(), alternative_gene_name: "AT2G34580|T31E10.8|T31E10_8".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20061019".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:2062356".to_string() },
            /* 43 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:2062356".to_string(), db_object_symbol: "AT2G34580".to_string(), invert: "".to_string(), go_term: "GO:0005739".to_string(), reference: "TAIR:AnalysisReference:501780126".to_string(), evidence_code: "ISM".to_string(), additional_evidence: "".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "".to_string(), alternative_gene_name: "AT2G34580|AT2G34580.2".to_string(), gene_product_type: "protein".to_string(), taxon: "taxon:3702".to_string(), date: "20180831".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:gene:4515101221".to_string() },

            // AT4G30872
            /* 44 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(), go_term: "GO:0005575".to_string(), reference: "TAIR:Communication:1345790".to_string(), evidence_code: "ND".to_string(), additional_evidence: "".to_string(), aspect: Aspect::CellularComponent,  unique_gene_name: "AT4G30872".to_string(), alternative_gene_name: "AT4G30872".to_string(), gene_product_type: "RNA".to_string(), taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string(), },
            /* 45 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(), go_term: "GO:0003674".to_string(), reference: "TAIR:Communication:1345790".to_string(), evidence_code: "ND".to_string(), additional_evidence: "".to_string(), aspect: Aspect::MolecularFunction,  unique_gene_name: "AT4G30872".to_string(), alternative_gene_name: "AT4G30872".to_string(), gene_product_type: "RNA".to_string(), taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string(), },
            /* 46 */ AnnotationRecord { db: "TAIR".to_string(), database_id: "locus:4515103469".to_string(), db_object_symbol: "AT4G30872".to_string(), invert: "".to_string(), go_term: "GO:0008150".to_string(), reference: "TAIR:Communication:1345790".to_string(), evidence_code: "ND".to_string(), additional_evidence: "".to_string(), aspect: Aspect::BiologicalProcess,  unique_gene_name: "AT4G30872".to_string(), alternative_gene_name: "AT4G30872".to_string(), gene_product_type: "RNA".to_string(), taxon: "taxon:3702".to_string(), date: "20090508".to_string(), assigned_by: "TAIR".to_string(), annotation_extension: "".to_string(), gene_product_form_id: "TAIR:locus:4515103469".to_string(), },
        ];

        static ref EVIDENCE_CODES: &'static [&'static str] = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
        static ref TEST_ANNOTATIONS: Vec<Annotation<'static>> = TEST_ANNOTATION_RECORDS.iter()
            .map(|record| Annotation::from_record(record, &EVIDENCE_CODES))
            .collect();
    }

    #[test]
    fn test_query_all() {
        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);
        let result = Query::All.execute(&index);

        // All of the genes from the input should appear in the query result
        assert!(TEST_GENES.iter().all(|gene| result.queried_genes.contains(gene)));

        // All of the annotations from the input should appear in the query result
        assert!(TEST_ANNOTATIONS.iter().all(|anno| result.queried_annotations.contains(anno)));
    }

    #[test]
    fn test_query_segment_bp_exp() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);
        let segment = Segment { aspect: BiologicalProcess, annotation_status: KnownExperimental };
        let result = segment.query(&index);

        let expected_genes_vec = vec![
            &TEST_GENES[0],
            &TEST_GENES[1],
            &TEST_GENES[2],
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &result.queried_genes);

        let expected_annotations_vec = vec![
            // AT5G48870
            &TEST_ANNOTATIONS[7],
            &TEST_ANNOTATIONS[9],

            // AT1G07060
            &TEST_ANNOTATIONS[14],
            &TEST_ANNOTATIONS[17],

            // AT4G34200
            &TEST_ANNOTATIONS[24],
            &TEST_ANNOTATIONS[25],
            &TEST_ANNOTATIONS[34],
            &TEST_ANNOTATIONS[39],
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &result.queried_annotations);
    }

    #[test]
    fn test_query_segment_mf_other() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);
        let segment = Segment { aspect: MolecularFunction, annotation_status: KnownOther };
        let result = segment.query(&index);

        let expected_genes_vec = vec![
            &TEST_GENES[0],
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &result.queried_genes);
        let expected_annotations_vec = vec![
            // AT5G48870
            &TEST_ANNOTATIONS[8],
            &TEST_ANNOTATIONS[10],
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &result.queried_annotations);
    }

    #[test]
    fn test_query_union() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);

        let segment_a = Segment { aspect: BiologicalProcess, annotation_status: KnownExperimental };
        let segment_b = Segment { aspect: MolecularFunction, annotation_status: KnownOther };
        let segment_c = Segment { aspect: CellularComponent, annotation_status: KnownOther };
        let query = Query::Union(vec![segment_a, segment_b, segment_c]);
        let results = query.execute(&index);

        let expected_genes_vec = vec![
            &TEST_GENES[0],
            &TEST_GENES[1],
            &TEST_GENES[2],
            &TEST_GENES[3],
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations_vec = vec![
            // AT5G48870
            &TEST_ANNOTATIONS[7],
            &TEST_ANNOTATIONS[8],
            &TEST_ANNOTATIONS[9],
            &TEST_ANNOTATIONS[10],

            // AT1G07060
            &TEST_ANNOTATIONS[12],
            &TEST_ANNOTATIONS[14],
            &TEST_ANNOTATIONS[15],
            &TEST_ANNOTATIONS[16],
            &TEST_ANNOTATIONS[17],
            &TEST_ANNOTATIONS[18],

            // AT4G34200
            &TEST_ANNOTATIONS[24],
            &TEST_ANNOTATIONS[25],
            &TEST_ANNOTATIONS[34],
            &TEST_ANNOTATIONS[39],

            // AT2G34580
            &TEST_ANNOTATIONS[40],
            &TEST_ANNOTATIONS[43],
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &results.queried_annotations);
    }

    #[test]
    fn test_query_union_unknowns() {
        use {Aspect::*, AnnotationStatus::*};

        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);
        let segment_a = Segment { aspect: BiologicalProcess, annotation_status: Unknown };
        let segment_b = Segment { aspect: MolecularFunction, annotation_status: Unknown };
        let segment_c = Segment { aspect: CellularComponent, annotation_status: Unknown };
        let query = Query::Union(vec![segment_a, segment_b, segment_c]);
        let results = query.execute(&index);

        let expected_genes_vec = vec![
            &TEST_GENES[3],
            &TEST_GENES[4],
        ];
        let expected_genes: HashSet<_> = expected_genes_vec.into_iter().collect();
        assert_eq!(&expected_genes, &results.queried_genes);

        let expected_annotations_vec = vec![
            // AT2G34580
            &TEST_ANNOTATIONS[41],
            &TEST_ANNOTATIONS[42],

            // AT4G30872
            &TEST_ANNOTATIONS[44],
            &TEST_ANNOTATIONS[45],
            &TEST_ANNOTATIONS[46],
        ];
        let expected_annotations: HashSet<_> = expected_annotations_vec.into_iter().collect();
        assert_eq!(&expected_annotations, &results.queried_annotations);
    }

    #[test]
    fn test_query_all_is_ordered() {
        let index = Index::new(&*TEST_GENES, &*TEST_ANNOTATIONS);
        let query = Query::All;
        let results = query.execute(&index);

        // Test that annotations are in the same order
        results.annotations_iter().zip(TEST_ANNOTATIONS.iter())
            .for_each(|(actual, expected)| assert_eq!(actual, expected));
    }
}
