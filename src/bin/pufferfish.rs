use clap::{Parser, Subcommand, ValueEnum};
use env_logger;
use lib::{
    fasta::{FastaReader, ParseError, Record},
    pufferfish,
    pufferfish::DefaultPufferfishIndex,
};
use postcard::{from_bytes, to_stdvec};
use std::{
    fs::File,
    io,
    io::{Read, Write},
};

#[derive(Parser, Debug)]
#[command(name = "pufferfish", version)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(ValueEnum, Clone, Debug, Default)]
pub enum Strategy {
    #[default]
    Default,
    Better,
}

#[derive(Debug, Subcommand)]
enum Command {
    Index {
        #[arg(long, short = 'f', value_delimiter = ' ', num_args = 1..)]
        files: Vec<String>,

        #[arg(long, short = 'k')]
        k: usize,

        #[arg(long, short = 'o')]
        out_file: String,

        #[arg(value_enum, default_value_t)]
        strategy: Strategy,
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

fn main() -> io::Result<()> {
    env_logger::init();

    let args = Args::parse();

    match &args.command {
        Command::Index {
            files,
            k,
            out_file,
            strategy,
        } => {
            let mut out_file = File::create(out_file).expect("Coud not create output file");

            let strategy = match strategy {
                Strategy::Default => pufferfish::Strategy::Default,
                Strategy::Better => pufferfish::Strategy::Better,
            };

            let records: Vec<(String, String)> = files
                .into_iter()
                .map(|file| {
                    let in_file = File::open(file).expect("File not found");
                    let mut reader = FastaReader::new(in_file);
                    reader
                        .records()
                        .map(|record| match record {
                            Ok(Record {
                                identifier,
                                sequence,
                            }) => (identifier, sequence),
                            Err(ParseError::IoError(err)) => panic!("Parse error: {err:?}"),
                            Err(ParseError::FormatError(err)) => panic!("Format error: {err:?}"),
                        })
                        .collect::<Vec<_>>()
                })
                .flatten()
                .collect();

            let index = DefaultPufferfishIndex::new(*k, records, strategy);

            let bytes: Vec<u8> = to_stdvec(&index).unwrap();
            out_file.write_all(&bytes).expect("Failed to write bytes");
        }
        Command::Query {
            index,
            query_file,
            out_file,
        } => {
            let mut index_file = File::open(index).expect("File not found");
            let query_file = File::open(query_file).expect("File not found");
            let mut bytes: Vec<u8> = Vec::new();
            index_file.read_to_end(&mut bytes).expect("Read failed");

            let mut out_file = File::create(out_file).expect("Coud not create output file");

            let index: DefaultPufferfishIndex = from_bytes(&bytes).unwrap();

            let mut reader = FastaReader::new(query_file);
            for record in reader.records() {
                let Record {
                    identifier,
                    sequence,
                } = record.expect("failed to read record");
                write!(out_file, "{identifier}\n")?;
                let found_colors = index.query(sequence);
                write!(out_file, "found in: {found_colors:?}\n")?;
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

    Ok(())
}
