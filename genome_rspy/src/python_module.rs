use pyo3::prelude::*;
use pyo3::PyErr;
use pyo3::exceptions::{PyIOError, PyValueError};


mod ex {
    use pyo3::{create_exception, exceptions::PyValueError};

    create_exception!(genome_rspy, HeatDeathError, PyValueError);
    create_exception!(genome_rspy, IndexSizeError, PyValueError);
    create_exception!(genome_rspy, InvalidNucleotideCharacterError, PyValueError);
}

#[pymodule]
mod genome_rspy {
    // See User's guide: https://pyo3.rs/v0.22.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.22.1/exception.html

    use std::path::Path;
    use crate::{GenomeError, OwnedHitRecord};

    use super::*;

    #[pymodule_export]
    use ex::{HeatDeathError, IndexSizeError, InvalidNucleotideCharacterError};

    /// Record of a read-chromosome match
    #[pyclass]
    #[allow(dead_code)]
    struct HitRecord {
        /// name of the read that was found
        #[pyo3(get, set)]
        read_name: String,

        /// where in the chromosome the match starts (it ends at chromosome_start + kmer_length)
        #[pyo3(get, set)]
        chromosome_start: usize
    }

    impl From<GenomeError> for PyErr {
        fn from(err: GenomeError) -> PyErr {
            match err {
                GenomeError::IOError(_) => PyIOError::new_err(err.to_string()),
                GenomeError::ThreadPoolError(_) => PyValueError::new_err(err.to_string()),
                GenomeError::HeatDeath(_,_) => HeatDeathError::new_err(err.to_string()),
                GenomeError::IndexSize(_) => IndexSizeError::new_err(err.to_string()),
                GenomeError::InvalidNucleotideCharacter => InvalidNucleotideCharacterError::new_err(err.to_string()),
            }
        }
    }

    impl From<OwnedHitRecord> for HitRecord {
        fn from(record: OwnedHitRecord) -> HitRecord {
            HitRecord {
                read_name: record.read_name,
                chromosome_start: record.chromosome_start
            }
        }
    }

    /// Given a path to a fasta file and a directory to store chromosomes in, loads and indexes chromosomes and returns a list of the names of the loaded chromosomes
    #[pyfunction]
    fn parse_chromosomes(
        py: Python<'_>,
        fasta_path: &str,
        storage_dir: &str,
        index_size: usize,
        num_threads: usize
    ) -> PyResult<Vec<(String, Option<PyErr>)>> {
        py.allow_threads(|| {
            crate::create_pool(num_threads)?.install(|| {
                let out = crate::parse_chromosomes(&Path::new(fasta_path), &Path::new(storage_dir), index_size)?;

                let out = out.into_iter()
                    .map(|(name, maybe_err)| (name, maybe_err.map(|v| v.into())))
                    .collect();

                Ok(out)
            })
        })
    }

    /// Search a pre-indexed chromosome for `kmer_length`-long subsequences of the sequences held in the fasta file specified by `reads_path`
    #[pyfunction]
    fn search_chromosome(
        py: Python<'_>,
        reads_path: &str,
        chromosome_name: &str,
        storage_dir: &str,
        kmer_length: usize,
        num_threads: usize
    ) -> PyResult<Vec<HitRecord>> {
        py.allow_threads(|| {
            crate::create_pool(num_threads)?.install(|| {
                let out = crate::load_mutants_and_search_chromosome(
                    &Path::new(reads_path),
                    chromosome_name,
                    &Path::new(storage_dir),
                    kmer_length
                )?;

                let out: Vec<HitRecord> = out.into_iter()
                    .map(|v| v.into())
                    .collect();
                
                Ok(out)
            })
        })
    }
}
