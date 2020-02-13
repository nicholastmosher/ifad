use clap::{App, Arg, ArgMatches};
use std::io::BufReader;
use ifad::MetadataReader;

fn app<'a, 'b>() -> clap::App<'a, 'b> {
    App::new("ifad")
        .arg(Arg::with_name("genes")
            .long("--genes")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("annotations")
            .long("--annotations")
            .required(true)
            .takes_value(true))
}

fn main() {
    let matches = app().get_matches();

    match run(&matches) {
        Ok(()) => (),
        Err(e) => eprintln!("{:?}", e),
    }
}

fn run(args: &ArgMatches) -> Result<(), String> {
    let genes_path = args.value_of("genes").expect("should have genes arg");
    let annos_path = args.value_of("annotations").expect("should have annotations arg");

    let mut genes_file = std::fs::File::open(genes_path)
        .map_err(|e| format!("failed to open genes file: {:?}", e))?;
    let mut gene_reader = MetadataReader::new(BufReader::new(&mut genes_file));
    let gene_records = ifad::GeneRecord::parse_from(&mut gene_reader)
        .map_err(|e| format!("failed to parse gene records: {:?}",e ))?;
    let gene_metadata = gene_reader.metadata().expect("should capture gene metadata");
    println!("Got {} genes", gene_records.len());
    println!("Gene metadata:\n{}", gene_metadata);

    let mut annos_file = std::fs::File::open(annos_path)
        .map_err(|e| format!("failed to open annotations file: {:?}", e))?;
    let mut anno_reader = MetadataReader::new(BufReader::new(&mut annos_file));
    let anno_records = ifad::AnnotationRecord::parse_from(&mut anno_reader)
        .map_err(|e| format!("failed to parse annotation records: {:?}", e))?;
    let anno_metadata = anno_reader.metadata().expect("should capture annotation metadata");
    println!();
    println!("Got {} annotations", anno_records.len());
    println!("Annotation metadata:\n{}", anno_metadata);

    Ok(())
}
