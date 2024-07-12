use std::fs;
use std::{io, path::Path};
use std::fmt::Display;

use packed_genome::{DeSerializable, PackedSequence, SimplePackedSequence, StandardIndexedPackedSequence};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelBridge, ParallelIterator};
use tqdm::tqdm;

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct Fasta(Vec<(String, String)>);

impl Fasta {
    // Fasta parsing code from reddit u/hch12908 (https://www.reddit.com/r/rust/comments/r5je0y/comment/hmnchrj/)
    pub fn parse_fasta(input: &str) -> Result<Fasta, ()> {
        // This iterator produces a Result<(name, sequence)> chain
        let sequences = input.split('>')
            // the first element will be empty (eg ">James\nAAAA>John\nGGGG" splits
            // into ["", "James\nAAAA", "John\nGGGG"])
            .skip(1)
            .par_bridge()
            .map(|section| {
                let mut lines = section.lines();
                let identifier = lines.next();
                let dna = lines.collect::<String>();

                match identifier {
                    Some(id) => Ok((id.into(), dna)),
                    None => Err(())
                }
            }).collect::<Vec<_>>();

        let mut fasta = Vec::new();

        for seq in sequences {
            fasta.push(seq?);
        }

        Ok(Fasta(fasta))
    }

    pub fn write_fasta(&self) -> String {
        // Create an empty String to store the FASTA formatted output
        let mut fasta_output = String::new();

        // Iterate over each (header, sequence) pair in the vector
        for (header, sequence) in &self.0 {
            // Append the header, starting with '>'
            fasta_output.push_str(&format!(">{}\n", header));
            // Append the sequence
            fasta_output.push_str(&format!("{}\n", sequence));
        }

        // Return the FASTA formatted string
        fasta_output
    }

    /// Clean data, removing leading/trailing 'N' characters in sequences and converting any
    /// 'N'-containing sequences into their first 'N'-free sequence of length >= kmer_length,
    /// dropping any sequences that don't fit these requirements
    pub fn clean(&mut self, kmer_length: usize) {
        let took = std::mem::take(&mut self.0);
        self.0 = tqdm(took.into_iter())
            .desc(Some("cleaning fasta"))
            .map(|(label, seq)| (label, seq.trim_matches('N').to_owned()))
            .map(|(label, seq)| {
                let seq = seq.split('N')
                    .filter(|v| v.len() >= kmer_length)
                    .map(|v| v.to_owned())
                    .next();
                return (label, seq);
            })
            .filter_map(|(label, seq)| Some((label, seq?)))
            .collect();
    }
}

#[derive(Clone)]
pub struct PackedFasta(Vec<(String, SimplePackedSequence)>);

impl From<PackedFasta> for Fasta {
    fn from(value: PackedFasta) -> Self {
        Fasta(
            value.0.into_iter()
            .map(
                |(label, seq)| (label, seq.str())
            )
            .collect()
        )
    }
}
impl From<Fasta> for PackedFasta {
    fn from(value: Fasta) -> Self {
        PackedFasta(
            value.0.into_iter()
            .map(
                |(label, seq)| (label, SimplePackedSequence::new(&seq))
            )
            .collect()
        )
    }
}

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Chromosome {
    name: String,
    raw_sequence: String,
    leading_ns: usize,
    seq: StandardIndexedPackedSequence
}

impl Chromosome {
    fn new(name: &str, raw_sequence: &str) -> Result<Chromosome, String> {
        let raw_sequence = raw_sequence.to_ascii_uppercase();
        let raw_sequence = raw_sequence.trim_end_matches('N');
        let pre_len = raw_sequence.len();

        let raw_sequence = raw_sequence.trim_start_matches('N');
        let leading_ns = pre_len - raw_sequence.len();

        let raw_sequence = raw_sequence.replace('N', "A");

        if raw_sequence.chars().any(|c| c != 'A' && c != 'C' && c != 'T' && c != 'G') {
            return Err("Invalid character, must be all A, C, T, G, or N".to_owned());
        }

        let seq = match StandardIndexedPackedSequence::new(&raw_sequence, 8) {
            Some(s) => s,
            None => return Err("Failed to create an IndexedPackedSequence with chunk length 8".to_owned())
        };

        let name = name.to_owned();

        Ok(Chromosome { name, raw_sequence, leading_ns, seq })
    }
}

impl Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Chromosome({} leading Ns {} sequence length {})", self.name, self.leading_ns, self.seq.len())
    }
}

macro_rules! iores {
    ($e:expr) => {
        match $e {
            Ok(v) => Ok(v),
            Err(e) => Err(io::Error::new(io::ErrorKind::Other, format!("{:?}", e)))
        }
    };
}

/// Given a path to a fasta file and a directory to store chromosomes in, returns a list of the
/// names of loaded chromosomes
pub fn parse_chromosomes(fasta_path: &Path, storage_dir: &Path) -> io::Result<Vec<String>> {
    fs::create_dir_all(storage_dir)?;

    println!("Reading file");
    let fasta_contents = fs::read_to_string(fasta_path)?;

    println!("Parsing fasta");
    let parsed = iores!(Fasta::parse_fasta(&fasta_contents))?;
    println!("Parsing chromosomes");

    return Ok(parsed.0.into_par_iter()
        .filter(|(label, _)| label == "5 dna:chromosome chromosome:GRCh38:5:1:181538259:1 REF")
        .map(|(label, seq)| {
            println!("> Parsing chromosome '{}'", &label);
            let chromosome = Chromosome::new(&label, &seq);
            drop(seq);

            match chromosome {
                Ok(c) => {
                    println!("> Writing chromosome '{}'", &label);
                    let file_name = format!("{:x}.chr", md5::compute(&label));
                    let file_path = storage_dir.join(file_name);
                    println!("> Serializing and compressing '{}'", &label);
                    let compressed = c.serialize_and_compress();

                    drop(c);

                    let write_result = match compressed {
                        Ok(comp) => {
                            println!("> Writing '{}' to disk", &label);
                            fs::write(file_path, comp)
                        },
                        Err(e) => Err(e)
                    };

                    match write_result {
                        Ok(()) => println!("> Saved chromosome '{}' successfully", &label),
                        Err(e) => println!("> Error: Failed to save chromosome '{}': {}", &label, e)
                    };
                },
                Err(e) => {
                    println!("> Error: Failed to load chromosome '{}': {}", label, e);
                }
            }

            return label.to_owned();
        })
        .collect()
    );
}

pub fn load_mutants(mutants_path: &Path) -> io::Result<PackedFasta> {
    println!("> Reading mutants");
    let contents = fs::read_to_string(mutants_path)?;

    println!("> Parsing mutants");
    let mut mutants = iores!(Fasta::parse_fasta(&contents))?;
    mutants.clean(31);

    packed_genome::disable_tqdm();
    let out = mutants.into();
    packed_genome::enable_tqdm();

    return Ok(out);
}

pub struct HitRecord {
    mutant_name: String,
    chromosome_name: String,
    // NOTE: the index written here is
    // (chromosome.raw_sequence.idx(mutant) + chromosome.leading_ns)
    chromosome_start: usize,
    chromosome_end: usize
}

pub fn search_chromosome(chromosome_name: &str, storage_dir: &Path, mutants: PackedFasta) -> io::Result<Vec<HitRecord>> {
    let file_name = format!("{:x}.chr", md5::compute(&chromosome_name));
    let file_path = storage_dir.join(file_name);

    println!("> Reading chromosome '{}'", &chromosome_name);
    let contents = fs::read(file_path)?;

    println!("> Decompressing and deserializing chromosome '{}'", &chromosome_name);
    let chromosome = iores!(Chromosome::decompress_and_deserialize(&contents))?;

    println!("> Searching chromosome '{}'", &chromosome_name);

    let out = tqdm(mutants.0.iter())
        .desc(Some("searching"))
        .map(move |(label, seq)| {
            seq.subsections(31)
            .map(move |subs| (label, subs))
        })
        .flatten()
        //.collect::<Vec<_>>()
        //.into_iter()
        .par_bridge()
        .map(|(label, subseq)| {
            let mut occurrences: Vec<HitRecord> = Vec::new();
            let mut start = 0;

            while start < chromosome.seq.len() - 31 {
                start  = match chromosome.seq.find_bounded(&subseq, Some(start), None) {
                    Some(found) => {
                        occurrences.push(HitRecord {
                            mutant_name: label.clone(),
                            chromosome_name: chromosome_name.to_owned(),
                            chromosome_start: found + chromosome.leading_ns,
                            chromosome_end: found + chromosome.leading_ns + subseq.len()
                        });
                        found + 1
                    },
                    None => break
                };
            }

            return occurrences;
        })
        .flatten_iter()
        .collect();

    return Ok(out);
}

#[cfg(test)]
mod tests;
