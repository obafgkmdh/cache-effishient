use std::{
    io::{BufRead, BufReader, Read},
    iter::Iterator,
    slice,
};

#[derive(Debug)]
pub struct Record {
    pub identifier: String,
    pub sequence: String,
}

#[derive(Debug)]
pub enum ParseError {
    IoError(std::io::Error),
    FormatError(String),
}

impl From<std::io::Error> for ParseError {
    fn from(err: std::io::Error) -> Self {
        ParseError::IoError(err)
    }
}

impl From<String> for ParseError {
    fn from(err: String) -> Self {
        ParseError::FormatError(err)
    }
}

pub struct FastaReader<R: Read> {
    reader: BufReader<R>,
}

pub struct Records<'a, R: Read> {
    reader: &'a mut FastaReader<R>,
}

impl<R: Read> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
        }
    }
}

impl<R: Read> FastaReader<R> {
    pub fn next_record(&mut self) -> Result<Option<Record>, ParseError> {
        let mut header_byte: u8 = 0;

        if self
            .reader
            .read_exact(slice::from_mut(&mut header_byte))
            .is_err()
        {
            // Reached EOF
            return Ok(None);
        }

        // Header should start with '>'
        (header_byte == b'>')
            .then_some(())
            .ok_or_else(|| format!("Expected header byte b'>', got {:?}", header_byte as char))?;

        let mut identifier = String::new();
        self.reader.read_line(&mut identifier)?;
        // Trim line ending
        if identifier.ends_with('\n') {
            identifier.pop();
            if identifier.ends_with('\r') {
                identifier.pop();
            }
        }

        let mut sequence = String::new();
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.first().is_none_or(|b| *b == b'>') {
                // Reached header line or EOF
                break;
            }
            self.reader.read_line(&mut sequence)?;
            // Trim line ending
            if sequence.ends_with('\n') {
                sequence.pop();
                if sequence.ends_with('\r') {
                    sequence.pop();
                }
            }
        }

        Ok(Some(Record {
            identifier,
            sequence,
        }))
    }

    pub fn records(&mut self) -> Records<'_, R> {
        Records { reader: self }
    }
}

impl<R: Read> Iterator for Records<'_, R> {
    type Item = Result<Record, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next_record().transpose()
    }
}
