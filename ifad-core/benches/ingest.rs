use criterion::{Criterion, criterion_group, criterion_main};
use ifad::{MetadataReader, AnnotationRecord};
use std::io::{Cursor, Read, BufReader};

const ANNOTATIONS_1000: &'static str = include_str!("truncated_1_000_tair.gaf");

fn run_metadata_annotations(data: &str) {
    let mut reader = MetadataReader::new(Cursor::new(data));
    let mut data_output = String::new();
    reader.read_to_string(&mut data_output).unwrap();
    let _metadata = reader.metadata().unwrap();
}

fn metadata_annotations_benchmark(c: &mut Criterion) {
    c.bench_function("metadata 1000 lines", |b| b.iter(|| run_metadata_annotations(ANNOTATIONS_1000)));
}

fn run_parse_annotations(data: &str) {
    let mut data_reader = Cursor::new(data);
    let _records = AnnotationRecord::parse_from(&mut data_reader).unwrap();
}

fn parse_annotations_benchmark(c: &mut Criterion) {
    let mut reader_1000 = MetadataReader::new(BufReader::new(Cursor::new(ANNOTATIONS_1000)));
    let mut data_1000 = String::with_capacity(ANNOTATIONS_1000.len());
    reader_1000.read_to_string(&mut data_1000).unwrap();

    c.bench_function("parse annotations 1000 lines", |b| b.iter(|| run_parse_annotations(&data_1000)));
}

criterion_group!(benches,
    metadata_annotations_benchmark,
    parse_annotations_benchmark);
criterion_main!(benches);
