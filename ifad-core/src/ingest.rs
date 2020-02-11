use std::io::Read;
use serde::{Deserialize, Serialize};
use crate::Aspect;

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct AnnotationRecord {
    pub db: String,
    pub database_id: String,
    pub db_object_symbol: String,
    pub invert: String,
    pub go_term: String,
    pub reference: String,
    pub evidence_code: String,
    pub additional_evidence: String,
    pub aspect: Aspect,
    pub unique_gene_name: String,
    pub alternative_gene_name: String,
    pub gene_product_type: String,
    pub taxon: String,
    pub date: String,
    pub assigned_by: String,
    pub annotation_extension: String,
    pub gene_product_form_id: String,
}

impl AnnotationRecord {
    pub fn parse_from<R: Read>(reader: R) -> Result<Vec<Self>, ()> {
        let mut csv_reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .flexible(true)
            .from_reader(reader);

        let results: Vec<_> = csv_reader.deserialize::<Self>()
            .filter_map(Result::ok)
            .collect();

        Ok(results)
    }
}

#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct GeneRecord {
    pub gene_id: String,
    pub gene_product_type: String,
}

impl GeneRecord {
    pub fn parse_from<R: Read>(reader: R) -> Result<Vec<Self>, ()> {
        let mut csv_reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .flexible(true)
            .from_reader(reader);

        let results: Vec<_> = csv_reader.deserialize::<Self>()
            .filter_map(Result::ok)
            .collect();

        Ok(results)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_annotation() {
        let annotation_string = "TAIR	locus:2031476	ENO1		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT1G74030	AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10	protein	taxon:3702	20190907	InterPro		TAIR:locus:2031476";
        let mut reader = Cursor::new(annotation_string);
        let annotations = AnnotationRecord::parse_from(&mut reader).unwrap();
        let expected = AnnotationRecord {
            db: "TAIR".to_string(),
            database_id: "locus:2031476".to_string(),
            db_object_symbol: "ENO1".to_string(),
            invert: "".to_string(),
            go_term: "GO:0000015".to_string(),
            reference: "TAIR:AnalysisReference:501756966".to_string(),
            evidence_code: "IEA".to_string(),
            additional_evidence: "InterPro:IPR000941".to_string(),
            aspect: Aspect::CellularComponent,
            unique_gene_name: "AT1G74030".to_string(),
            alternative_gene_name: "AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10".to_string(),
            gene_product_type: "protein".to_string(),
            taxon: "taxon:3702".to_string(),
            date: "20190907".to_string(),
            assigned_by: "InterPro".to_string(),
            annotation_extension: "".to_string(),
            gene_product_form_id: "TAIR:locus:2031476".to_string(),
        };
        assert_eq!(vec![expected], annotations);
    }

    #[test]
    fn test_parse_gene() {
        let gene_row = "AT1G01010	protein_coding";
        let mut reader = Cursor::new(gene_row);
        let genes = GeneRecord::parse_from(&mut reader).unwrap();
        let expected = GeneRecord {
            gene_id: "AT1G01010".to_string(),
            gene_product_type: "protein_coding".to_string(),
        };
        assert_eq!(vec![expected], genes);
    }
}
