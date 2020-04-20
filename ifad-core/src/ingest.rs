use std::io::{Read, BufRead, Cursor, Error};
use serde::{Deserialize, Serialize};
use crate::Aspect;

pub struct MetadataReader<B> {
    reader: B,
    metadata: String,
    header: String,
    metadata_finished: bool,
    buffer: Cursor<String>,
}

impl<B: BufRead> MetadataReader<B> {
    pub fn new(reader: B) -> MetadataReader<B> {
        MetadataReader {
            reader,
            metadata: String::new(),
            header: String::new(),
            metadata_finished: false,
            buffer: Cursor::new(String::new()),
        }
    }

    pub fn metadata(&self) -> Option<&str> {
        if !self.metadata_finished { return None; }
        Some(&self.metadata)
    }

    pub fn header(&self) -> Option<&str> {
        if !self.metadata_finished { return None; }
        Some(&self.header)
    }
}

impl<B: BufRead> Read for MetadataReader<B> {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, Error> {
        let len = self.buffer.read(buf)?;
        if len != 0 { return Ok(len); }

        // If we finished reading from the metadata section, just
        // continue reading from the inner reader.
        if self.metadata_finished {
            return self.reader.read(buf);
        }

        loop {

            // Read a line into the internal buffer
            let len = self.reader.read_line(&mut self.buffer.get_mut())?;
            if len == 0 { return Ok(0); }

            let line = self.buffer.get_ref().trim_start();
            if line.is_empty() {
                self.metadata.push('\n');
                self.buffer.get_mut().clear();
            } else if line.starts_with('!') {
                // Add this line to the metadata
                self.metadata.push_str(self.buffer.get_ref());
                self.buffer.get_mut().clear();
            } else {
                // Mark that the metadata section has ended
                self.metadata_finished = true;

                // This line must be the header
                self.buffer.read_to_string(&mut self.header)?;
                self.buffer.get_mut().clear();
                return self.reader.read(buf);
            }
        }
    }
}

#[cfg_attr(test, derive(Clone))]
#[derive(Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
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

        let mut records = Vec::new();
        let mut row = csv::StringRecord::new();
        while csv_reader.read_record(&mut row).unwrap() {
            let record: AnnotationRecord = row.deserialize(None).unwrap();
            records.push(record);
        }

        Ok(records)
    }
}

#[derive(Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(test, derive(Clone))]
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

        let mut records = Vec::new();
        let mut row = csv::StringRecord::new();
        while csv_reader.read_record(&mut row).unwrap() {
            let record: GeneRecord = row.deserialize(None).unwrap();
            records.push(record);
        }

        Ok(records)
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

    #[test]
    fn test_metadata_reader() {
        let input = r"


!gaf-version: 2.1
!


!Generated by GO Central


!
!Date Generated by GOC: 2019-10-07


DB	DB Object ID	DB Object Symbol	Qualifier	GO ID	DB:Reference (JDB:Reference)	Evidence Code	With (or) From	Aspect	DB Object Name	DB Object Type	Taxon	Date	Assigned By	Annotation Extension	Gene Product Form ID
TAIR	locus:2031476	ENO1		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT1G74030	AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10	protein	taxon:3702	20190907	InterPro		TAIR:locus:2031476
TAIR	locus:2043067	ENOC		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G29560	AT2G29560|ENOC|ENO3|cytosolic enolase|enolase 3|F16P2.6|F16P2_6	protein	taxon:3702	20190408	InterPro		TAIR:locus:2043067
TAIR	locus:2044851	LOS2		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G36530	AT2G36530|LOS2|ENO2|LOW EXPRESSION OF OSMOTICALLY RESPONSIVE GENES 2|enolase 2|F1O11.16|F1O11_16	protein	taxon:3702	20190408	InterPro		TAIR:locus:2044851
TAIR	locus:2032970	AT1G25260		GO:0000027	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR033867	P	AT1G25260	AT1G25260|F4F7.35|F4F7_35	protein	taxon:3702	20190404	InterPro		TAIR:locus:2032970";
        let input_cursor = Cursor::new(input);
        let mut metadata_reader = MetadataReader::new(input_cursor);
        let mut output = String::new();
        let _ = metadata_reader.read_to_string(&mut output);
        let metadata_output = metadata_reader.metadata().unwrap();
        let header_output = metadata_reader.header().unwrap();

        let expected_metadata = r"


!gaf-version: 2.1
!


!Generated by GO Central


!
!Date Generated by GOC: 2019-10-07


";
        assert_eq!(metadata_output, expected_metadata);

        let expected_header = "DB	DB Object ID	DB Object Symbol	Qualifier	GO ID	DB:Reference (JDB:Reference)	Evidence Code	With (or) From	Aspect	DB Object Name	DB Object Type	Taxon	Date	Assigned By	Annotation Extension	Gene Product Form ID\n";
        assert_eq!(header_output, expected_header);

        let expected_body = r"TAIR	locus:2031476	ENO1		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT1G74030	AT1G74030|ENO1|enolase 1|F2P9.10|F2P9_10	protein	taxon:3702	20190907	InterPro		TAIR:locus:2031476
TAIR	locus:2043067	ENOC		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G29560	AT2G29560|ENOC|ENO3|cytosolic enolase|enolase 3|F16P2.6|F16P2_6	protein	taxon:3702	20190408	InterPro		TAIR:locus:2043067
TAIR	locus:2044851	LOS2		GO:0000015	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR000941	C	AT2G36530	AT2G36530|LOS2|ENO2|LOW EXPRESSION OF OSMOTICALLY RESPONSIVE GENES 2|enolase 2|F1O11.16|F1O11_16	protein	taxon:3702	20190408	InterPro		TAIR:locus:2044851
TAIR	locus:2032970	AT1G25260		GO:0000027	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR033867	P	AT1G25260	AT1G25260|F4F7.35|F4F7_35	protein	taxon:3702	20190404	InterPro		TAIR:locus:2032970";
        assert_eq!(output, expected_body);
    }
}
