use clap::{Parser, Subcommand};

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

        #[arg(long, short = 'o')]
        out_file: String,
    },
    Query {
        #[arg(long, short = 'i')]
        index: String,

        #[arg(long, short = 'o')]
        out_file: String,
    },
}

fn main() {
    let args = Args::parse();

    match &args.command {
        Command::Index { file, out_file } => {
            println!("got args: {:?} {:?}", file, out_file);
            todo!();
        },
        Command::Query { index, out_file } => {
            println!("got args: {:?} {:?}", index, out_file);
            todo!();
        },
    }
}
