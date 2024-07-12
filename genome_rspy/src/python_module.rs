use pyo3::prelude::*;
use pyo3::PyErr;
use pyo3::exceptions::{PyIOError, PyValueError};


mod ex {
    use pyo3::{create_exception, exceptions::PyException};

    create_exception!(genome_rspy, HeatDeathError, PyException);
}

#[pymodule]
mod genome_rspy {
    // See User's guide: https://pyo3.rs/v0.22.1/
    // mutable example (no return) see https://github.com/PyO3/rust-numpy
    // https://pyo3.rs/v0.22.1/exception.html

    use std::path::Path;
    use crate::GenomeError;

    use super::*;

    #[pymodule_export]
    use ex::HeatDeathError;

    impl From<GenomeError> for PyErr {
        fn from(err: GenomeError) -> PyErr {
            match err {
                GenomeError::IOError(_) => PyIOError::new_err(err.to_string()),
                GenomeError::ThreadPoolError(_) => PyValueError::new_err(err.to_string()),
                GenomeError::HeatDeath(_, _) => HeatDeathError::new_err(err.to_string()),
            }
        }
    }

    /// Given a path to a fasta file and a directory to store chromosomes in, loads and indexes chromosomes and returns a list of the names of the loaded chromosomes
    #[pyfunction]
    fn parse_chromosomes(
        py: Python<'_>,
        fasta_path: &str,
        storage_dir: &str,
        num_threads: usize
    ) -> PyResult<Vec<String>> {
        py.allow_threads(|| {
            crate::create_pool(num_threads)?.install(|| {
                Ok(crate::parse_chromosomes(&Path::new(fasta_path), &Path::new(storage_dir))?)
            })
        })
    }
}
