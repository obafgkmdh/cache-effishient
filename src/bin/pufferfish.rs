use clap::{Parser, Subcommand};
use lib::{
    fasta::{FastaReader, ParseError, Record},
    pufferfish::DefaultPufferfishIndex,
};
use postcard::{from_bytes, to_stdvec};
use std::{
    fs::File,
    io::{Read, Write},
};

#[derive(Parser, Debug)]
#[command(name = "naive", version)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Index {
        #[arg(long, short = 'f')]
        file: String,

        #[arg(long, short = 'k')]
        k: usize,

        #[arg(long, short = 'o')]
        out_file: String,
    },
    Query {
        #[arg(long, short = 'i')]
        index: String,

        #[arg(long, short = 'q')]
        query_file: String,

        #[arg(long, short = 'o', default_value = "/dev/stdout")]
        out_file: String,
    },
    Inspect {
        #[arg(long, short = 'i')]
        index: String,
    },
}

fn main() {
    let args = Args::parse();

    match &args.command {
        Command::Index { file, k, out_file } => {
            let in_file = File::open(file).expect("File not found");
            let mut out_file = File::create(out_file).expect("Coud not create output file");
            let mut reader = FastaReader::new(in_file);
            let sequences: Vec<String> = reader
                .records()
                .map(|record| match record {
                    Ok(Record {
                        identifier: _,
                        sequence,
                    }) => sequence,
                    Err(ParseError::IoError(err)) => panic!("Parse error: {err:?}"),
                    Err(ParseError::FormatError(err)) => panic!("Format error: {err:?}"),
                })
                .collect();

            let index = DefaultPufferfishIndex::new(*k, sequences);

            let bytes: Vec<u8> = to_stdvec(&index).unwrap();
            out_file.write_all(&bytes).expect("Failed to write bytes");
        }
        Command::Query {
            index,
            query_file,
            out_file: _,
        } => {
            let mut index_file = File::open(index).expect("File not found");
            let query_file = File::open(query_file).expect("File not found");
            let mut bytes: Vec<u8> = Vec::new();
            index_file.read_to_end(&mut bytes).expect("Read failed");

            let index: DefaultPufferfishIndex = from_bytes(&bytes).unwrap();

            let mut reader = FastaReader::new(query_file);
            for record in reader.records() {
                let Record {
                    identifier,
                    sequence,
                } = record.expect("failed to read record");
                let found = index.query(sequence);
                println!("{identifier}: {found}");
            }
        }
        Command::Inspect { index } => {
            let mut index_file = File::open(index).expect("File not found");
            let mut bytes: Vec<u8> = Vec::new();
            index_file.read_to_end(&mut bytes).expect("Read failed");

            let index: DefaultPufferfishIndex = from_bytes(&bytes).unwrap();

            index.print_stats();
        }
    }
}
