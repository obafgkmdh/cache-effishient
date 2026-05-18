// Generates a random series of queries based upon specified flags

use clap::{Parser, Subcommand};
use env_logger;
use log::debug;
use lib::fasta::{FastaReader, ParseError, Record};
use std::path::Path;
use std::{collections::BTreeSet, fs::File, io::Write};

use rand::distr::{Distribution, slice::Choose};
use rand::prelude::*;

#[derive(Parser, Debug)]
#[command(name = "genq", version)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    Random {
        #[arg(long, short = 'm')]
        min: usize,

        #[arg(long, short = 'M')]
        max: usize,

        #[arg(long, short = 'n')]
        num: usize,

        #[arg(long, short = 'o')]
        out_file: String,
    },
    FromGenomes {
        #[arg(long, short = 'f', value_delimiter = ' ', num_args = 1..)]
        genome_files: Vec<String>,

        #[arg(long, short = 'm')]
        min: usize,

        #[arg(long, short = 'M')]
        max: usize,

        #[arg(long, short = 'n')]
        num: usize,

        #[arg(long, short = 'p')]
        perfile_min: usize,

        #[arg(long, short = 'o')]
        out_file: String,
    },
    FromGenes {
        #[arg(long, short = 'f', value_delimiter = ' ', num_args = 1..)]
        gene_files: Vec<String>,

        #[arg(long, short = 'm')]
        min: usize,

        #[arg(long, short = 'M')]
        max: usize,

        #[arg(long, short = 'n')]
        num: usize,

        #[arg(long, short = 'o')]
        out_file: String,
    },
}

fn main() {
    env_logger::init();

    let mut rng = rand::rng();
    let args = Args::parse();

    match args.command {
        Command::Random {
            min,
            max,
            num,
            out_file,
        } => {
            assert!(min <= max, "invalid bounds ({}, {})", min, max);
            assert!(num > 0, "need num > 0 (got {})", num);
            let mut out_file = File::create(out_file).expect("Coud not create output file");
            let alphabet = ['A', 'C', 'T', 'G'];
            let dist = Choose::new(&alphabet).unwrap();

            for i in 1..=num {
                let sample_size = rng.random_range(min..=max);
                let string: String = dist.sample_iter(&mut rng).take(sample_size).collect();

                write!(out_file, ">id:{i} size:{sample_size}\n{string}\n")
                    .expect("output write failed");
            }
        }
        Command::FromGenomes {
            genome_files,
            min,
            max,
            num,
            perfile_min,
            out_file,
        } => {
            assert!(min <= max, "invalid bounds ({}, {})", min, max);
            assert!(num > 0, "need num > 0 (got {})", num);

            let n_files = genome_files.len();
            assert!(
                n_files * perfile_min <= num,
                "Cannot sample from {} files with given minimum {perfile_min} and number requested {num}!",
                n_files
            );

            let mut out_file = File::create(out_file).expect("Coud not create output file");
            let remaining = num - n_files * perfile_min;

            // determine allocation of remaining samples to each file
            // we do this by inserting n_files - 1 "dividers" among the remaining samples
            let mut random_positions: BTreeSet<usize> = BTreeSet::new();
            let n_slots = remaining + n_files - 1;
            while random_positions.len() < n_files - 1 {
                random_positions.insert(rng.random_range(0..n_slots));
            }
            let mut pick_from_file: Vec<usize> = Vec::with_capacity(n_files);
            let mut last = 0;
            for position in random_positions {
                pick_from_file.push(perfile_min + position - last);
                last = position + 1;
            }
            pick_from_file.push(perfile_min + n_slots - last);

            debug!("count from each file: {:?}", pick_from_file);

            for (file_path, mut n_samples) in
                genome_files.into_iter().zip(pick_from_file.into_iter())
            {
                let file = File::open(&file_path).expect("File not found");
                let file_name = Path::new(&file_path).file_name().unwrap();

                let mut reader = FastaReader::new(file);
                let records: Vec<Record> = reader
                    .records()
                    .map(|record| match record {
                        Ok(r) => r,
                        Err(ParseError::IoError(err)) => panic!("Parse error: {err:?}"),
                        Err(ParseError::FormatError(err)) => panic!("Format error: {err:?}"),
                    })
                    .collect();

                while n_samples > 0 {
                    let Record {
                        identifier,
                        sequence,
                    } = records.choose(&mut rng).unwrap();
                    let sample_size = rng.random_range(min..=max);

                    if sample_size > sequence.len() {
                        continue;
                    }

                    let start = rng.random_range(0..=(sequence.len() - sample_size));
                    let sample = &sequence[start..(start + sample_size)];
                    write!(
                        out_file,
                        ">file: {file_name:?} size:{sample_size} record: {identifier}\n{sample}\n"
                    )
                    .expect("output write failed");

                    n_samples -= 1;
                }
            }
        }
        Command::FromGenes {
            gene_files,
            min,
            max,
            num,
            out_file,
        } => {
            assert!(min <= max, "invalid bounds ({}, {})", min, max);
            assert!(num > 0, "need num > 0 (got {})", num);
            let mut out_file = File::create(out_file).expect("Coud not create output file");

            let mut files: Vec<File> = Vec::new();

            for file_name in gene_files {
                let file = File::open(file_name).expect("File not found");
                eprintln!("{:?}", file);
                files.push(file);
            }
        }
    }
}
