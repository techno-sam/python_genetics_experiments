use std::fs;
use std::{io, path::Path};
use std::fmt::Display;

use packed_genome::{PackedSequence, SimplePackedSequence, StandardIndexedPackedSequence};

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
            .map(|section| {
                let mut lines = section.lines();
                let identifier = lines.next();
                let dna = lines.collect::<String>();

                match identifier {
                    Some(id) => Ok((id.into(), dna)),
                    None => Err(())
                }
            });

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
        self.0 = took.into_iter()
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
    fn new(name: &str, raw_sequence: &str) -> Chromosome {
        let raw_sequence = raw_sequence.trim_end_matches('N');
        let pre_len = raw_sequence.len();

        let raw_sequence = raw_sequence.trim_start_matches('N');
        let leading_ns = pre_len - raw_sequence.len();

        let raw_sequence = raw_sequence.replace('N', "A");

        let seq = StandardIndexedPackedSequence::new(&raw_sequence, 8)
            .expect("Given chunk length should be valid");

        let name = name.to_owned();

        Chromosome { name, raw_sequence, leading_ns, seq }
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

    let fasta_contents = fs::read_to_string(fasta_path)?;

    let parsed = iores!(Fasta::parse_fasta(&fasta_contents))?;

    todo!();
}

#[cfg(test)]
mod tests;
