use crate::{Aspect, AnnotationStatus};
use crate::ingest::AnnotationRecord;

#[derive(Debug, Hash, Eq, PartialEq)]
pub struct Annotation {
    pub db: String,
    pub database_id: String,
    pub db_object_symbol: String,
    pub invert: String,
    pub go_term: String,
    pub reference: String,
    pub evidence_code: String,
    pub additional_evidence: String,
    pub aspect: Aspect,
    pub annotation_status: AnnotationStatus,
    pub gene_names: Vec<String>,
    pub gene_product_type: String,
    pub taxon: String,
    pub date: String,
    pub assigned_by: String,
    pub annotation_extension: String,
    pub gene_product_form_id: String,
}

impl From<AnnotationRecord> for Annotation {
    fn from(record: AnnotationRecord) -> Self {
        let gene_names = vec![record.unique_gene_name]
            .into_iter()
            .chain(record.alternative_gene_name
                .split("|")
                .map(ToOwned::to_owned))
            .collect();

        Annotation {
            db: record.db,
            database_id: record.database_id,
            db_object_symbol: record.db_object_symbol,
            invert: record.invert,
            go_term: record.go_term,
            reference: record.reference,
            evidence_code: record.evidence_code,
            additional_evidence: record.additional_evidence,
            aspect: record.aspect,
            annotation_status: unimplemented!(),
            gene_names,
            gene_product_type: record.gene_product_type,
            taxon: record.taxon,
            date: record.date,
            assigned_by: record.assigned_by,
            annotation_extension: record.annotation_extension,
            gene_product_form_id: record.gene_product_form_id,
        }
    }
}

#[derive(Debug)]
pub struct Gene {
    pub gene_id: String,
    pub gene_product_type: String,
}
