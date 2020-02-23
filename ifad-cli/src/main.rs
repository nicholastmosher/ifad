use clap::{App, Arg, ArgMatches, Values, AppSettings};
use std::io::BufReader;
use ifad::{MetadataReader, Annotation, Gene, Index, Segment, GafExporter, Query};
use std::convert::TryFrom;

fn app<'a, 'b>() -> clap::App<'a, 'b> {

    let segment_validator = |segment: String| {
        let splits: Vec<&str> = segment.split(',').collect();
        if splits.len() != 2 { return Err("segments must be written as ASPECT,STATUS".to_string()); }
        if let [aspect, status] = splits[..] {
            if !&["F", "C", "P"].contains(&aspect) {
                return Err("aspect must be one of F, C, or P".to_string());
            }
            if !&["EXP", "OTHER", "UNKNOWN", "UNANNOTATED"].contains(&status) {
                return Err("status must be one of EXP, OTHER, UNKNOWN, or UNANNOTATED".to_string());
            }
            return Ok(());
        }
        unreachable!()
    };

    App::new("ifad")
        .setting(AppSettings::DeriveDisplayOrder)
        .arg(Arg::with_name("genes")
            .help("the file to read genes from (e.g. gene-types.txt")
            .long("--genes")
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations")
            .help("the file to read annotations from (e.g. tair.gaf)")
            .long("--annotations")
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("genes_out")
            .help("the file to write queried genes to (e.g. gene-types_F-EXP.txt")
            .long("--genes-out")
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations_out")
            .help("the file to write queried annotations to (e.g. tair_F-EXP.gaf)")
            .long("--annotations-out")
            .require_equals(true)
            .takes_value(true))
        .arg(Arg::with_name("query")
            .help("the type of query")
            .long("--query")
            .possible_values(&["union", "intersection"])
            .default_value("union")
            .require_equals(true))
        .arg(Arg::with_name("segment")
            .help("a segment to use in the query, given as ASPECT,STATUS (e.g. F,EXP or C,OTHER)")
            .multiple(true)
            .long("--segment")
            .require_equals(true)
            .takes_value(true)
            .validator(segment_validator))
}

fn main() {
    let matches = app().get_matches();

    match run(&matches) {
        Ok(()) => (),
        Err(e) => eprintln!("{}", e),
    }
}

struct Config<'a> {
    genes_path: &'a str,
    annos_path: &'a str,
    genes_out: &'a str,
    annos_out: &'a str,
    query: &'a str,
    segments: Values<'a>,
}

impl Config<'_> {
    fn from_args<'a>(args: &'a ArgMatches) -> Option<Config<'a>> {
        let genes_path = args.value_of("genes")?;
        let annos_path = args.value_of("annotations")?;
        let genes_out = args.value_of("genes_out")?;
        let annos_out = args.value_of("annotations_out")?;
        let query = args.value_of("query")?;
        let segments = args.values_of("segment")?;
        Some(Config { genes_path, annos_path, genes_out, annos_out, query, segments })
    }
}

fn run(args: &ArgMatches) -> Result<(), String> {
    let maybe_config = Config::from_args(args);
    let config = match maybe_config {
        Some(options) => options,
        None => {
            app().print_help().unwrap();
            return Ok(());
        },
    };

    let segments: Vec<Segment> = config.segments.map(|segment| {
        let split: Vec<&str> = segment.split(',').collect();
        let segment = (split[0], split[1]);
        Segment::try_from(segment).expect("should convert segment arg to Segment")
    }).collect();

    let mut genes_file = std::fs::File::open(config.genes_path)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut gene_reader = MetadataReader::new(BufReader::new(&mut genes_file));
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}",e ))?;
    let gene_metadata = gene_reader.metadata().expect("should capture gene metadata");
    let gene_headers = gene_reader.header().expect("should get gene headers");

    let mut annos_file = std::fs::File::open(config.annos_path)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata");
    let anno_headers = anno_reader.header().expect("should capture annotation header");

    let genes: Vec<Gene> = gene_records.iter()
        .map(|record| Gene::from_record(record))
        .collect();

    let experimental_evidence = &["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP"];
    let annotations: Vec<Annotation> = anno_records.iter()
        .map(|record| Annotation::from_record(record, experimental_evidence))
        .collect();

    let index: Index = Index::new(&genes, &annotations);
    let query = match config.query {
        "union" => Query::Union(segments),
        // "intersection" => Query::Intersection(segments),
        "intersection" => return Err("Intersection queries are not yet implemented!".to_string()),
        _ => unreachable!(),
    };

    eprintln!("Executing query: {:?}", query);
    let result = query.execute(&index);

    let mut genes_out = std::fs::File::create(config.genes_out)
        .map_err(|e| format!("failed to create genes output file: {:?}", e))?;
    let mut genes_exporter = GafExporter::new(
        gene_metadata.to_string(),
        gene_headers.to_string(),
        result.genes_iter().map(|gene| gene.record));
    genes_exporter.write_all(&mut genes_out).expect("should write genes file");

    let mut annotations_out = std::fs::File::create(config.annos_out)
        .map_err(|e| format!("failed to create annotations output file: {:?}", e))?;
    let mut annotations_exporter = GafExporter::new(
        anno_metadata.to_string(),
        anno_headers.to_string(),
        result.annotations_iter().map(|anno| anno.record));
    annotations_exporter.write_all(&mut annotations_out)
        .map_err(|e| format!("failed to export data as GAF: {:?}", e))?;

    Ok(())
}
