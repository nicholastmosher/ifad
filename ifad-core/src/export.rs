use std::io::Write;
use serde::Serialize;

pub struct GafExporter<I: Iterator> {
    metadata: String,
    header: String,
    record_iter: I,
}

impl<T, I: Iterator<Item=T>> GafExporter<I>
    where T: Serialize
{
    pub fn new(
        metadata: String,
        header: String,
        record_iter: I,
    ) -> GafExporter<I> {
        GafExporter { metadata, header, record_iter }
    }

    pub fn write_all<W: Write>(&mut self, mut writer: W) -> std::io::Result<()> {
        write!(&mut writer, "{}", self.metadata)?;
        write!(&mut writer, "{}", self.header)?;

        let mut csv_writer = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(writer);
        for record in &mut self.record_iter {
            csv_writer.serialize(record)?;
        }
        csv_writer.flush()?;
        Ok(())
    }
}
