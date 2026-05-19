# Cache-Effishient

To fetch the e-coli genomes used for testing, run the following script:
```sh
#!/bin/sh
wget https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.190/Mitarbeiter/beller/CPM2015/ecoli.tar.gz
mkdir -p data/ecoli
tar -xf ecoli.tar.gz -C data/ecoli/
cd data/ecoli/
rm 6.*.fasta 7.*.fasta 12.*.fasta 17.*.fasta 19.*.fasta 25.*.fasta 26.*.fasta 38.*.fasta 46.*.fasta 52.*.fasta 54.*.fasta 56.*.fasta 58.*.fasta
```

For index creation, run the following:
```sh
cargo run --release --bin pufferfish -- index -f data/ecoli/*.fasta -k 23 --strategy greedy -o data/ecoli.bin
```
One can also set the `RUST_LOG` environment variable to `debug` or `trace` for verbose messages. For Pufferfish behavior, use `--strategy default` instead.

To generate queries, run the following:
```sh
cargo run --release --bin genq -- from-genomes -f data/ecoli/*.fasta -m 1000 -M 5000 -n 100000 -p 400 -o data/randgenomes.fa
```

To run queries against the index, run the following:
```sh
cargo run --release --bin pufferfish -- query -i data/ecoli.bin -q data/randgenomes.fa -H 2
```
The `-H` flag controls the optimization level. `-H 0` does no optimization, `-H 1` emulates the behavior of Pufferfish, and `-H 2` does the most optimization.
