use std::path::Path;

use genome_rspy::{load_mutants, parse_chromosomes, search_chromosome, search_chromosome_general};
use packed_genome::{indexed_packed_sequence, IndexedPackedSequence, PackedSequence, SimplePackedSequence};

fn main() {
    let seq = &"TTCAGCCGAACTGGATCAGAACGGATAGTATATCGACACCGCACTCCGGCTTCATTCATGCACG
        ATACCACGATTTGAACTAGAAGAAGCTCGAAGCATTCCCTCTGGTTATAACTCATATTCAAATTCGAATTGA
        AAAGGAAATTATACGGATAAGTGACGCGTGTTAGTTACATCAGTGCAGCTCGTGAATTTTGCCAGTACACGG
        TCTCGGATCCCATGGTGCCGCAGCGCCGGAGGCCGTGACGTGCAGTTCTACTCATATCTCAGCGCTTGAAAA
        CTTCCAGATAGCACAAGGCTACACGGTGCACCGTCAGGTATTTTTATGACGACAACCATCTGTCAAGTTGCA
        TGGCCGCTGCATATCTCTTGCAGGCTCTACAATGAGGGGGACCGACCAGCTGACATTTCGGAGTTGCGAGTC
        CCGTGCTCTAAAACCTTTGCTTCGCGCCAGGTAAGACAGCTGAGCGCATGACTGTGCGTGTCCGCTTAGAGT
        ATCTCCACTGCGATATTAGACACCACTGTTACCCGCGTACTGCGGCCGACCCAATACATACACAGAAAACTA
        CCCACAGAGCTCACTCTGCGCGTGGGAATTAGCTGATCAAGTTCACGATAGTCGTACTACTCCATTCCGGTG
        GCACGGGGAGATGTCGTTAATGGACTTCTCCGCAATACTAACGCTAATTTTGCCTTTGAGTATTCATAAAAG
        ACCCGACGTGTACATCGAGCAATCATAGGTTCAGATGAATCCGTCGCAACCGTTAAATTTGTGGGGTTGGGC
        ACCGTTTGGATATCCGATCGGGTTAAGCTGTCGCTCCGGACGCTGCCGATAGACACGCGCAAATGTATATTA
        GAGCGTATTGCACTTGTTGGCGGAGACGTGCAGGTTGTAACGCGTTATGCAGACGGGGTGGCTCCCCGCTCT
        ATTAAAGTGGTTCTCGTAATGGTCTTCGCAGACCACGGTGCCTGCGGCCTAATCGATTTTTTCACTGCCTTG"
        .replace("\n", "").replace(" ", "");

    println!("Sequence: {seq}\n");

    let packed = indexed_packed_sequence!(seq, 4);
    println!("Packed: {}", packed);

    let mut start = None;

    loop {
        start = packed.find_bounded_str("GGGGACCG", start, None);
        
        start = match start {
            Some(s) => {
                println!("Found needle starting at {}", s);
                Some(s+1)
            },
            None => break
        };
    }

    for v in SimplePackedSequence::new("AGTTAAGCTGTCGCTCCGGACGCTGCCGATAGACACGCGCAAATGTATATTA").subsections(7) {
        println!("Subsection: {}", v);
    }

    for v in SimplePackedSequence::new("ACT").subsections(7) {
        println!("Subsection: {}", v);
    }

    let parse: bool = false;

    if parse {
        let r = parse_chromosomes(
            Path::new("/home/sam/PycharmProjects/PythonADClassWorkspace/python_genetics_experiments/scratch/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
            Path::new("/home/sam/PycharmProjects/PythonADClassWorkspace/python_genetics_experiments/scratch/rust_checkpoints")
        );

        match r {
            Ok(v) => {
                println!("Parsed {} chromosomes:", v.len());
                for name in v {
                    println!("\t{}", name);
                }
            },
            Err(e) => println!("Failed to parse chromosomes: {}", e)
        };
    } else {
        let mutants = load_mutants(
            Path::new("/home/sam/PycharmProjects/PythonADClassWorkspace/python_genetics_experiments/scratch/gget_mutate_out_mutant_reference.fa")
        ).expect("Failed to load mutants");

        let r = search_chromosome_general(
            "5 dna:chromosome chromosome:GRCh38:5:1:181538259:1 REF",
            Path::new("/home/sam/PycharmProjects/PythonADClassWorkspace/python_genetics_experiments/scratch/rust_checkpoints"),
            &mutants,
            31
        );

        match r {
            Ok(v) => {
                println!("Searched chromosome, got {} hits", v.len());
                for hit in v {
                    println!("\t{:?}", hit);
                }
            },
            Err(e) => println!("Failed to search chromosome: {}", e)
        };
    }
}
