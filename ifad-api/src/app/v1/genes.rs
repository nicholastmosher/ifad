use arc_swap::ArcSwap;
use futures::StreamExt;
use serde::{Deserialize, Serialize};
use actix_web::error::PayloadError;
use actix_web::{HttpResponse, web};

use ifad::{Query, Segment, StreamingGafExporter};
use crate::app::AppData;

#[derive(Debug, Serialize, Deserialize)]
pub enum Filter {
    #[serde(rename = "all")]
    All,
    #[serde(rename = "include_protein")]
    ProteinCoding,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum Strategy {
    #[serde(rename = "union")]
    Union,
    #[serde(rename = "intersection")]
    Intersection,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum Format {
    #[serde(rename = "json")]
    Json,
    #[serde(rename = "gene-csv")]
    GeneCSV,
    #[serde(rename = "gaf")]
    Gaf,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GenesQuery {
    filter: Filter,
    strategy: Strategy,
    format: Format,
}

#[derive(Debug, Serialize, Deserialize)]
struct CountResponse<'a> {
    gene_count: usize,
    annotation_count: usize,
    gene_metadata: &'a str,
    annotation_metadata: &'a str,
}

pub async fn read(
    state: web::Data<ArcSwap<AppData>>,
    query: web::Query<GenesQuery>,
    json: web::Json<Vec<Segment>>,
) -> Result<HttpResponse, ()> {
    let appdata = state.load_full();
    let index = appdata.index.clone();

    let query = query.into_inner();
    let segments = json.into_inner();
    let query_object = match query.strategy {
        Strategy::Union => Query::Union(segments),
        Strategy::Intersection => Query::Intersection(segments),
    };
    let query_result = query_object.execute(index);

    match query.format {
        Format::Json => {
            let data = CountResponse {
                gene_count: query_result.gene_count(),
                annotation_count: query_result.annotation_count(),
                gene_metadata: &appdata.gene_metadata,
                annotation_metadata: &appdata.anno_metadata,
            };
            Ok(HttpResponse::Ok().json(data))
        }
        Format::Gaf => {
            let stream = StreamingGafExporter::new(
                appdata.anno_metadata.to_string(),
                appdata.anno_headers.to_string(),
                query_result.iter_annotations().map(|anno| anno.record)
            ).map(|result| result.map_err(|_| PayloadError::EncodingCorrupted));
            Ok(HttpResponse::Ok().streaming(stream))
        }
        Format::GeneCSV => {
            let stream = StreamingGafExporter::new(
                appdata.gene_metadata.to_string(),
                appdata.gene_headers.to_string(),
                query_result.iter_genes().map(|gene| gene.record)
            ).map(|result| result.map_err(|_| PayloadError::EncodingCorrupted));
            Ok(HttpResponse::Ok().streaming(stream))
        }
    }
}
