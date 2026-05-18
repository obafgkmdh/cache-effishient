// Generates a random series of queries based upon specified flags

use clap::{Parser, Subcommand};
use env_logger;
use lib::{
    fasta::{FastaReader, ParseError, Record},
    pufferfish::DefaultPufferfishIndex,
};
use postcard::{from_bytes, to_stdvec};
use std::{
    fs::File,
    io::{Read, Write},
};

use rand::Rng;
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

struct Alphabet;

impl Distribution<char> for Alphabet {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> char {
        *b"ATCG".choose(rng).unwrap() as char
    }
}

fn random_string(n: usize) -> String {
    let mut rng = rand::rng();
    (0..n).map(|_| rng.sample(Alphabet)).collect()
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
            out_file,
        } => {
            let mut files: Vec<File> = Vec::new();

            for file_name in genome_files {
                let file = File::open(file_name).expect("File not found");
                println!("{:?}", file);
                files.push(file);
            }
        }
        Command::FromGenes {
            gene_files,
            min,
            max,
            num,
            out_file,
        } => {
            let mut files: Vec<File> = Vec::new();

            for file_name in gene_files {
                let file = File::open(file_name).expect("File not found");
                println!("{:?}", file);
                files.push(file);
            }
        }
    }
}
