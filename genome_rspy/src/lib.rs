#[cfg(test)]
mod tests;

mod python_module;

use std::fs;
use std::{io, path::Path};
use std::fmt::Display;

use packed_genome::{DeSerializable, PackedSequence, SimplePackedSequence, StandardIndexedPackedSequence};
use rayon::iter::{IntoParallelIterator, ParallelBridge, ParallelIterator};
use rayon::ThreadPoolBuildError;
use thiserror::Error;
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
            value.0.into_par_iter()
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
    fn new(name: &str, raw_sequence: &str, index_size: usize) -> Result<Chromosome, GenomeError> {
        let raw_sequence = raw_sequence.to_ascii_uppercase();
        let raw_sequence = raw_sequence.trim_end_matches('N');
        let pre_len = raw_sequence.len();

        let raw_sequence = raw_sequence.trim_start_matches('N');
        let leading_ns = pre_len - raw_sequence.len();

        let raw_sequence = raw_sequence.replace('N', "A");

        if raw_sequence.chars().any(|c| c != 'A' && c != 'C' && c != 'T' && c != 'G') {
            return Err(GenomeError::InvalidNucleotideCharacter);
        }

        let seq = match StandardIndexedPackedSequence::new(&raw_sequence, index_size) {
            Some(s) => s,
            None => return Err(GenomeError::IndexSize(index_size))
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
pub fn parse_chromosomes(fasta_path: &Path, storage_dir: &Path, index_size: usize) -> Result<Vec<(String, Option<GenomeError>)>, GenomeError> {
    fs::create_dir_all(storage_dir)?;

    println!("Reading file");
    let fasta_contents = fs::read_to_string(fasta_path)?;

    println!("Parsing fasta");
    let parsed = iores!(Fasta::parse_fasta(&fasta_contents))?;
    println!("Parsing chromosomes");

    return Ok(parsed.0.into_par_iter()
        .filter(|(label, _)| label == "5 dna:chromosome chromosome:GRCh38:5:1:181538259:1 REF") // FIXME:
                                                                                                // remove
                                                                                                // before
                                                                                                // release
        .map(|(label, seq)| {
            println!("> Parsing chromosome '{}'", &label);
            let chromosome = Chromosome::new(&label, &seq, index_size);
            drop(seq);

            let err = match chromosome {
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
                        Ok(()) => {
                            println!("> Saved chromosome '{}' successfully", &label);
                            None
                        },
                        Err(e) => {
                            println!("> Error: Failed to save chromosome '{}': {}", &label, e);
                            Some(e.into())
                        },
                    }
                },
                Err(e) => {
                    println!("> Error: Failed to load chromosome '{}': {}", label, e);
                    Some(e)
                }
            };

            return (label.to_owned(), err);
        })
        .collect()
    );
}

#[derive(Error, Debug)]
pub enum GenomeError {
    #[error("Cannot use index length {0} for {1}-mers, that will take until the heat death of the universe.")]
    /// (index length, kmer length)
    HeatDeath(usize, usize),

    #[error("Index size {0} is not a valid size. Valid sizes are 1, 2, 3, 4, 5, 6, 7, 8, 16, and 32.")]
    IndexSize(usize),

    #[error("Invalid nucleotide character, sequence must be made purely of A, C, T, G, or N characters.")]
    InvalidNucleotideCharacter,

    #[error(transparent)]
    IOError(#[from] std::io::Error),

    #[error(transparent)]
    ThreadPoolError(#[from] ThreadPoolBuildError)
}

pub fn create_pool(num_threads: usize) -> Result<rayon::ThreadPool, GenomeError> {
    return Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()?)
}

pub fn load_mutants_and_search_chromosome(
    mutants_path: &Path,
    chromosome_name: &str,
    storage_dir: &Path,
    kmer_length: usize
) -> Result<Vec<OwnedHitRecord>, GenomeError> {
    let mutants = load_mutants(mutants_path)?;

    let found = search_chromosome_general(chromosome_name, storage_dir, &mutants, kmer_length)?;

    Ok(found.into_iter().map(|v| v.into()).collect())
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

#[derive(Debug)]
pub struct HitRecord<'a> {
    read_name: &'a str,
    // NOTE: the index written here is
    // (chromosome.raw_sequence.idx(mutant) + chromosome.leading_ns)
    chromosome_start: usize
}

#[derive(Debug)]
pub struct OwnedHitRecord {
    read_name: String,
    chromosome_start: usize
}

impl<'a> From<HitRecord<'a>> for OwnedHitRecord {
    fn from(value: HitRecord<'a>) -> Self {
        OwnedHitRecord {
            read_name: value.read_name.to_owned(),
            chromosome_start: value.chromosome_start
        }
    }
}

pub fn search_chromosome_general<'a>(
    chromosome_name: &str,
    storage_dir: &Path,
    mutants: &'a PackedFasta,
    kmer_length: usize
) -> Result<Vec<HitRecord<'a>>, GenomeError> {
    if kmer_length == 31 {
        return search_chromosome(chromosome_name, storage_dir, mutants);
    }

    println!("Warning: searching with kmer length other than 31 is significantly slower");

    // This is exactly the same as below, but less generalized
    let file_name = format!("{:x}.chr", md5::compute(&chromosome_name));
    let file_path = storage_dir.join(file_name);

    println!("> Reading chromosome '{}'", &chromosome_name);
    let contents = fs::read(file_path)?;

    println!("> Deserializing chromosome '{}'", &chromosome_name);
    let chromosome = Chromosome::decompress_and_deserialize(&contents)?;

    println!("> Searching chromosome '{}'", &chromosome_name);

    let chunk_length = chromosome.seq.chunk_length();
    if chunk_length > kmer_length {
        return Err(GenomeError::HeatDeath(chunk_length, kmer_length));
    }

    let out = tqdm(mutants.0.iter())
        .desc(Some("searching"))
        .map(move |(label, seq)| {
            seq.subsections(kmer_length)
                .map(move |subs| (label, subs))
        })
        .flatten()
        .par_bridge()
        .map(|(label, subseq)| {
            let mut occurrences: Vec<HitRecord> = Vec::new();

            let all = chromosome.seq.find_all(&subseq);

            if let Some(v) = all {
                 for found in v {
                     occurrences.push(HitRecord {
                         read_name: &label,
                         chromosome_start: found + chromosome.leading_ns
                     });
                 }
            }

            return occurrences;
        })
        .flatten_iter()
        .take_any(5) // FIXME: remove this after done testing
        .collect();

    return Ok(out);
}

// TODO: Python API combined with `load_mutants` above
pub fn search_chromosome<'a>(
    chromosome_name: &str,
    storage_dir: &Path,
    mutants: &'a PackedFasta
) -> Result<Vec<HitRecord<'a>>, GenomeError> {
    let file_name = format!("{:x}.chr", md5::compute(&chromosome_name));
    let file_path = storage_dir.join(file_name);

    println!("> Reading chromosome '{}'", &chromosome_name);
    let contents = fs::read(file_path)?;

    println!("> Deserializing chromosome '{}'", &chromosome_name);
    let chromosome = Chromosome::decompress_and_deserialize(&contents)?;

    println!("> Searching chromosome '{}'", &chromosome_name);

    let chunk_length = chromosome.seq.chunk_length();
    if chunk_length > 31 {
        return Err(GenomeError::HeatDeath(chunk_length, 31));
    }

    let out = tqdm(mutants.0.iter())
        .desc(Some("searching"))
        .map(move |(label, seq)| {
            seq.subsections(31) // WARN: if this gets changed from 31, the line marked below MUST
                                // be changed as well. Perhaps you should use the general function?
            .map(move |subs| (label, subs))
        })
        .flatten()
        //.collect::<Vec<_>>()
        //.into_iter())
        .par_bridge()
        .map(|(label, subseq)| {
            let mut occurrences: Vec<HitRecord> = Vec::new();

            // WARN: if seq.subsections(31) gets changed to a different number, this call must be
            // changed to chromosome.seq.find_all(&subseq), which is much less optimized
            let all = unsafe {chromosome.seq.find_all_31mer(&subseq)};

            if let Some(v) = all {
                for found in v {
                    occurrences.push(HitRecord {
                        read_name: &label,
                        chromosome_start: found + chromosome.leading_ns,
                        //chromosome_end: found + chromosome.leading_ns + subseq.len()
                    });
                }
            }

            return occurrences;
        })
        .flatten_iter()
        .take_any(5) // FIXME: remove this after done testing
        .collect();

    return Ok(out);
}
