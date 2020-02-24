use futures::FutureExt;
use actix_web::{HttpServer, App, web, HttpResponse, Responder};
use ifad::{MetadataReader, Gene, Annotation, Index, GeneRecord, AnnotationRecord, Query, StreamingGafExporter};
use std::io::{BufReader, Read, BufRead};
use std::sync::Arc;
use arc_swap::ArcSwap;
use std::ops::Deref;
use actix_web::dev::HttpResponseBuilder;
use actix_web::http::StatusCode;
use futures::io::AllowStdIo;
use actix_web::error::PayloadError;

struct Config {
    genes_file: String,
    annotations_file: String,
}

impl Config {
    fn from_env() -> Option<Config> {
        let genes_file = std::env::var("GENES_FILE").ok()?;
        let annotations_file = std::env::var("ANNOTATIONS_FILE").ok()?;
        Some(Config { genes_file, annotations_file })
    }
}

fn main() {
    match run() {
        Ok(()) => return,
        Err(e) => eprintln!("{}", e),
    }
}

struct GeneData {
    gene_metadata: String,
    gene_headers: String,
    gene_records: Vec<GeneRecord>,
}

struct AnnotationData {
    anno_metadata: String,
    anno_headers: String,
    anno_records: Vec<AnnotationRecord>,
}

fn run() -> Result<(), String> {
    let config = Config::from_env()
        .ok_or("failed to read Config from environment")?;

    let mut genes_file = std::fs::File::open(config.genes_file)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut annos_file = std::fs::File::open(config.annotations_file)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;

    let mut swap = ArcSwap::new(Arc::new(None));
    swap_data(&mut swap, BufReader::new(&mut genes_file), BufReader::new(&annos_file))?;

    actix::System::new("ifad").block_on(server(swap).map(|_| ()));
    Ok(())
}

fn swap_data<GR, AR>(
    swap: &mut ArcSwap<Option<Index>>,
    mut genes: GR,
    mut annotations: AR,
) -> Result<(), String>
    where GR: BufRead, AR: BufRead,
{
    let mut gene_reader = MetadataReader::new(&mut genes);
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}", e))?;
    let gene_metadata = gene_reader.metadata().expect("should capture gene metadata").to_string();
    let gene_headers = gene_reader.header().expect("should get gene headers").to_string();
    let gene_data = Arc::new(GeneData { gene_metadata, gene_headers, gene_records });
    let genes: Vec<Gene> = gene_data.gene_records.into_iter()
        .map(|record| Gene::from_record(&record))
        .collect();
    let genes_arc = Arc::new(genes);

    let mut anno_reader = MetadataReader::new(&mut annotations);
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata").to_string();
    let anno_headers = anno_reader.header().expect("should capture annotation header").to_string();
    let anno_data = Arc::new(AnnotationData { anno_metadata, anno_headers, anno_records });
    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_data.anno_records.iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();
    let annos_arc = Arc::new(annotations);

    let index = Arc::new(Some(Index::new(&genes_arc, &annos_arc)));
    swap.store(index);
    Ok(())
}

async fn index(data: web::Data<ArcSwap<Option<Index<'static, 'static>>>>) -> impl Responder {
    use futures::StreamExt;

    let data: Option<&Index> = data.load().as_ref().as_ref();
    let index = match data {
        None => return HttpResponse::Ok().body("Data has not been loaded yet"),
        Some(index) => index,
    };

    let query_result = Query::All.execute(&index);
    let stream = StreamingGafExporter::new(
        "anno metadata".to_string(),
        "anno_header".to_string(),
        query_result.annotations_iter().map(|item| item.record))
        .map(|result| result.map_err(|_| PayloadError::EncodingCorrupted));

    HttpResponseBuilder::new(StatusCode::OK)
        .streaming(stream)
}

async fn server(data: ArcSwap<Option<Index<'static, 'static>>>) -> std::io::Result<()> {
    HttpServer::new(move || {
        App::new()
            .data(data)
            .route("/", web::get().to(index))
    })
    .bind("127.0.0.1:8000")?
    .run()
    .await
}
