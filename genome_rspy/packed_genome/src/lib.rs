use bincode;

//use bzip2::Compression;
//use bzip2::write::{BzEncoder, BzDecoder};

use ndarray::prelude::*;

use std::collections::HashMap;
use std::fmt::Display;
use std::io::{self/*, Write*/};
use std::ops::Deref;
use std::sync::atomic::{AtomicBool, Ordering};

static DO_TQDM: AtomicBool = AtomicBool::new(true);

pub fn enable_tqdm() {
    DO_TQDM.store(true, Ordering::SeqCst);
}

pub fn disable_tqdm() {
    DO_TQDM.store(false, Ordering::SeqCst);
}

#[cfg(feature = "progress")]
fn tqdm<'a, T: std::iter::Iterator + 'a>(v: T, desc: &str) -> Box<dyn std::iter::Iterator<Item = T::Item> + 'a> {
    use tqdm;

    if DO_TQDM.load(Ordering::SeqCst) {
        Box::new(tqdm::tqdm(v).desc(Some(desc)))
    } else {
        Box::new(v)
    }
}

#[cfg(not(feature = "progress"))]
fn tqdm<T: std::iter::Iterator>(v: T, desc: &str) -> T {
    v
}

fn pack(c: char) -> Option<u8> {
    match c {
        'A' => Some(0b00),
        'C' => Some(0b01),
        'T' => Some(0b10),
        'G' => Some(0b11),
        _ => None
    }
}

fn unpack(u: u8) -> Option<char> {
    match u {
        0b00 => Some('A'),
        0b01 => Some('C'),
        0b10 => Some('T'),
        0b11 => Some('G'),
        _ => None
    }
}

fn unpack_byte(u: u8) -> String {
    let a = (u & 0b11_00_00_00) >> 6;
    let b = (u & 0b00_11_00_00) >> 4;
    let c = (u & 0b00_00_11_00) >> 2;
    let d = u & 0b00_00_00_11;

    return [
        unpack(a).unwrap(),
        unpack(b).unwrap(),
        unpack(c).unwrap(),
        unpack(d).unwrap()
    ].into_iter().collect();
}

fn pack_string(seq: &str) -> Vec<u8> {
    return tqdm(seq.chars().collect::<Vec<_>>().chunks(4), "packing").map(|chunk| {
        let mut packed: u8 = 0;

        for j in 0..4 {
            packed <<= 2;
            if j < chunk.len() {
                let c = chunk[j];
                packed |= pack(c).unwrap();
            }
        }

        packed
    }).collect();
}

struct LocalString(String);
impl Deref for LocalString {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl From<String> for LocalString {
    fn from(value: String) -> Self {
        LocalString(value)
    }
}

trait FindBounded<Pat, T> where T: num::Integer {
    fn find_bounded(&self, pat: Pat, start: Option<T>, end: Option<T>) -> Option<T>;
}

pub trait DeSerializable: Sized + serde::Serialize + for <'a> serde::Deserialize<'a> {
    fn serialize_and_compress(&self) -> io::Result<Vec<u8>> {
        // Serialize the data to a vector of bytes
        let serialized_data = match bincode::serialize(self) {
            Ok(sd) => sd,
            Err(error) => return Err(io::Error::new(io::ErrorKind::Other, error))
        };

        // Compress the serialized data
        //let mut encoder = BzEncoder::new(Vec::new(), Compression::best());
        //encoder.write_all(&serialized_data)?;
        //let compressed_data = encoder.finish()?;
        let compressed_data = serialized_data;

        Ok(compressed_data)
    }

    fn decompress_and_deserialize(data: &[u8]) -> io::Result<Self> {
        // Decompress the data
        //let mut decoder = BzDecoder::new(Vec::new());
        //decoder.write_all(data)?;
        //let decompressed_data = decoder.finish()?;
        let decompressed_data = data;

        println!("> Deserializing {} bytes", decompressed_data.len());

        // Deserialize the decompressed data
        let deserialized_data: Self = match bincode::deserialize(&decompressed_data) {
            Ok(dd) => dd,
            Err(error) => return Err(io::Error::new(io::ErrorKind::Other, error))
        };

        Ok(deserialized_data)
    }
}

impl<T> DeSerializable for T where T: Sized + serde::Serialize + for <'a> serde::Deserialize<'a> {}

pub trait PackedSequence {
    fn new(seq: &str) -> Self where Self: Sized;

    fn len(&self) -> usize;

    fn get(&self, idx: usize) -> char;

    fn find(&self, pat: &impl PackedSequence) -> Option<usize> {
        self.find_bounded(pat, None, None)
    }
    fn find_str(&self, pat: &str) -> Option<usize> {
        self.find(&SimplePackedSequence::new(pat))
    }

    fn find_bounded(&self, pat: &impl PackedSequence, start: Option<usize>, end: Option<usize>) -> Option<usize>;
    fn find_bounded_str(&self, pat: &str, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        self.find_bounded(&SimplePackedSequence::new(pat), start, end)
    }

    fn find_all(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        let mut start = Some(0);

        let mut out = Vec::new();

        loop {
            start = match self.find_bounded(pat, start, None) {
                Some(s) => {
                    out.push(s);
                    Some(s + 1)
                },
                None => break
            };
        }

        return Some(out);
    }

    /// WARN: `pat` MUST be length 31 or this search WILL fail SILENTLY
    /// (overrides DO NOT have to check length)
    unsafe fn find_all_31mer(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        assert_eq!(31, pat.len());
        return self.find_all(pat);
    }

    fn contains(&self, other: &impl PackedSequence) -> bool {
        self.find(other).is_some()
    }
    fn contains_str(&self, other: &str) -> bool {
        self.find(&SimplePackedSequence::new(other)).is_some()
    }

    fn str(&self) -> String {
        (0..self.len()).map(|i| self.get(i)).collect()
    }

    fn subsections(&self, chunk_length: usize) -> impl Iterator<Item = SimplePackedSequence>;

    fn get_packed<'a>(&'a self) -> &'a Vec<u8>;

}

impl<T> From<&T> for LocalString where T: PackedSequence {
    fn from(value: &T) -> Self {
        value.str().into()
    }
}
/*impl<Sub, T> FindBounded<Sub, usize> for T where T: PackedSequence, Sub: PackedSequence {
    fn find_bounded(&self, pat: Sub, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        PackedSequence::find_bounded(self, pat, start, end)
    }
}*/

macro_rules! fbstr {
    ($type_name: ty) => {
        impl<'a> FindBounded<$type_name, usize> for String {
            fn find_bounded(&self, pat: $type_name, start: Option<usize>, end: Option<usize>) -> Option<usize> {
                let self_ = match end {
                    Some(end) => &self[..end],
                    None => &self
                };
                let (self_, offset) = match start {
                    Some(start) => (&self_[start..], start),
                    None => (self_, 0)
                };

                self_.find(pat).map(|v| v + offset)
            }
        }
    };
}
fbstr!(&str);
fbstr!(char);
fbstr!(&[char]);
fbstr!(&&str);
fbstr!(&String);

fn variants(v: Vec<u8>) -> [Vec<u8>; 4] {
    let mut variants = [v, Vec::new(), Vec::new(), Vec::new()];

    for i in 1..4 {
        let mut variant = Array::from_vec(variants[i - 1].clone());
        let mut copy = variant.clone();
        variant <<= 2;
        copy >>= 6;
        let mut slice = variant.slice_mut(s![..variant.len() -1]);
        slice |= &copy.slice(s![1..]);
        variants[i] = variant.to_vec();
    }

    variants
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct SimplePackedSequence {
    packed: Vec<u8>,
    len: usize
}

impl Display for SimplePackedSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SimplePackedSequence({} {})", self.len(), self.str())
    }
}

impl PackedSequence for SimplePackedSequence {
    fn new(seq: &str) -> SimplePackedSequence {
        let len = seq.len();
        let packed = pack_string(seq);

        return SimplePackedSequence { packed, len };
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.len
    }

    fn get(&self, idx: usize) -> char {
        let byte_idx = idx / 4;
        let shift = 6 - 2 * (idx % 4);
        return unpack((self.packed[byte_idx] >> shift) & 0b11).unwrap();
    }

    fn find_bounded(&self, sub: &impl PackedSequence, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        let self_ = self.str();
        let sub_ = sub.str();

        self_.find_bounded(&sub_, start, end)
    }

    fn subsections(&self, chunk_length: usize) -> impl Iterator<Item = SimplePackedSequence> {
        let variants = variants(self.packed.clone());

        let end = if chunk_length > self.len() {0} else {self.len() - chunk_length + 1};

        return (0..end).map(move |start| {
            let variant = &variants[start % 4];
            let packed = &variant[start/4 .. start/4 + (chunk_length+3)/4];

            SimplePackedSequence { packed: packed.into(), len: chunk_length }
        });
    }

    #[inline(always)]
    fn get_packed<'a>(&'a self) -> &'a Vec<u8> {
        &self.packed
    }
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct PreVariedPackedSequence {
    parent: SimplePackedSequence,
    variants: [Vec<u8>; 4]
}

impl Display for PreVariedPackedSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PreVariedPackedSequence({} {})", self.len(), self.str())
    }
}

impl PackedSequence for PreVariedPackedSequence {
    fn new(seq: &str) -> PreVariedPackedSequence {
        let parent = SimplePackedSequence::new(seq);

        {
            #[cfg(feature = "progress")]
            println!("> varying ...");
        }

        let variants = variants(parent.packed.clone());

        return PreVariedPackedSequence { parent, variants };
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.parent.len()
    }

    #[inline(always)]
    fn get(&self, idx: usize) -> char {
        self.parent.get(idx)
    }

    #[inline(always)]
    fn find_bounded(&self, pat: &impl PackedSequence, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        self.parent.find_bounded(pat, start, end)
    }

    fn subsections(&self, chunk_length: usize) -> impl Iterator<Item = SimplePackedSequence> {
        let end = if chunk_length > self.len() {0} else {self.len() - chunk_length + 1};
        return (0..end).map(move |start| {
            let variant = &self.variants[start % 4];
            let packed = &variant[start/4 .. start/4 + (chunk_length+3)/4];

            SimplePackedSequence { packed: packed.into(), len: chunk_length }
        });
    }

    #[inline(always)]
    fn get_packed<'a>(&'a self) -> &'a Vec<u8> {
        self.parent.get_packed()
    }
}

pub trait NucleotideKey<const C: usize> where Self: Sized + std::cmp::Eq + std::hash::Hash + Copy + Display {
    const BYTES_NEEDED: usize = (C + 3) / 4;
    const BYTES_NEEDED_FLOOR: usize = C / 4;

    fn str_to_key(seq: &str) -> Self {
        let seq = pack_string(seq);
        return Self::to_key(&seq);
    }

    fn to_key(seq: &[u8]) -> Self;
    
    unsafe fn to_key_unchecked(seq: &[u8]) -> Self {
        return Self::to_key(seq);
    }

    fn to_string(&self) -> String;
}

impl NucleotideKey<1> for u8 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 1);
        return seq[0] & 0b11_00_00_00;
    }

    fn to_string(&self) -> String {
        unpack_byte(*self)[..1].to_string()
    }
}
impl NucleotideKey<2> for u8 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 1);
        return seq[0] & 0b11_11_00_00;
    }

    fn to_string(&self) -> String {
        unpack_byte(*self)[..2].to_string()
    }
}
impl NucleotideKey<3> for u8 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 1);
        return seq[0] & 0b11_11_11_00;
    }

    fn to_string(&self) -> String {
        unpack_byte(*self)[..3].to_string()
    }
}
impl NucleotideKey<4> for u8 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 1);
        return seq[0];
    }

    fn to_string(&self) -> String {
        unpack_byte(*self)
    }
}

impl NucleotideKey<5> for u16 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 2);
        return ((seq[0] as u16) << 8) | ((seq[1] as u16) & 0b11_00_00_00);
    }

    fn to_string(&self) -> String {
        unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)[..1]
    }
}
impl NucleotideKey<6> for u16 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 2);
        return ((seq[0] as u16) << 8) | ((seq[1] as u16) & 0b11_11_00_00);
    }

    fn to_string(&self) -> String {
        unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)[..2]
    }
}
impl NucleotideKey<7> for u16 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 2);
        return ((seq[0] as u16) << 8) | ((seq[1] as u16) & 0b11_00_00_00);
    }

    fn to_string(&self) -> String {
        unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)[..3]
    }
}
impl NucleotideKey<8> for u16 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 2);
        return ((seq[0] as u16) << 8) | ((seq[1] as u16) & 0b11_00_00_00);
    }

    fn to_string(&self) -> String {
        unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)
    }
}

impl NucleotideKey<16> for u32 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 4);
        return ((seq[0] as u32) << 24) | ((seq[1] as u32) << 16) | ((seq[2] as u32) << 8) | (seq[3] as u32);
    }

    fn to_string(&self) -> String {
        unpack_byte((self >> 24) as u8) + &unpack_byte((self >> 16) as u8) + &unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)
    }
}
impl NucleotideKey<32> for u64 {
    fn to_key(seq: &[u8]) -> Self {
        assert!(seq.len() >= 8);
        return ((seq[0] as u64) << 56) | ((seq[1] as u64) << 48) | ((seq[2] as u64) << 40) | ((seq[3] as u64) << 32) | 
               ((seq[4] as u64) << 24) | ((seq[5] as u64) << 16) | ((seq[6] as u64) << 8) | (seq[7] as u64);
    }


    fn to_string(&self) -> String {
        unpack_byte((self >> 56) as u8) + &unpack_byte((self >> 48) as u8) + &unpack_byte((self >> 40) as u8) + &unpack_byte((self >> 32) as u8) +
       &unpack_byte((self >> 24) as u8) + &unpack_byte((self >> 16) as u8) + &unpack_byte((self >> 8) as u8) + &unpack_byte((self & 0xff) as u8)
    }
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
struct CheckedIndexedPackedSequence<K: NucleotideKey<C>, const C: usize> {
    parent: PreVariedPackedSequence,
    index: HashMap<K, Vec<usize>>,
    c: usize
}

impl<K: NucleotideKey<C>, const C: usize> TryFrom<CheckedIndexedPackedSequence<K, C>> for IndexedPackedSequence<K, C> {
    type Error = String;

    fn try_from(value: CheckedIndexedPackedSequence<K, C>) -> Result<Self, Self::Error> {
        if value.c == C {
            Ok(IndexedPackedSequence { parent: value.parent, index: value.index })
        } else {
            Err(format!("Check failed for IndexedPackedSequence, expected chunk length {}, got {}", C, value.c))
        }
    }
}

impl<K: NucleotideKey<C>, const C: usize> Into<CheckedIndexedPackedSequence<K, C>> for IndexedPackedSequence<K, C> {
    fn into(self) -> CheckedIndexedPackedSequence<K, C> {
        CheckedIndexedPackedSequence { parent: self.parent, index: self.index, c: C }
    }
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
#[serde(into = "CheckedIndexedPackedSequence<K, C>")]
#[serde(try_from = "CheckedIndexedPackedSequence<K, C>")]
pub struct IndexedPackedSequence<K: NucleotideKey<C>, const C: usize> {
    parent: PreVariedPackedSequence,
    index: HashMap<K, Vec<usize>>
}

impl<K: NucleotideKey<C>, const C: usize> Display for IndexedPackedSequence<K, C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "IndexedPackedSequence<{}>({} {} index length: {})", C, self.len(), self.str(), self.index.len())
    }
}

impl<K: NucleotideKey<C>, const C: usize> IndexedPackedSequence<K, C> {
    #[allow(dead_code)]
    fn debug_index(&self) {
        println!("Index:");
        for (k, v) in self.index.iter() {
            println!("\t{} {:?}", k.to_string(), v);
        }
    }

    #[allow(dead_code)]
    fn str_to_key(&self, seq: &str) -> K {
        return K::str_to_key(seq);
    }
}

impl<K: NucleotideKey<C>, const C: usize> PackedSequence for IndexedPackedSequence<K, C> {
    fn new(seq: &str) -> Self where Self: Sized {
        let parent = PreVariedPackedSequence::new(seq);
        let mut index: HashMap<K, Vec<usize>> = HashMap::new();

        for i in tqdm(0..(parent.len() - C + 1), "indexing") {
            let sample_from = &parent.variants[i % 4];
            let mut chunk = Vec::new();
            chunk.reserve_exact(K::BYTES_NEEDED);

            for j in 0..K::BYTES_NEEDED {
                chunk.push(sample_from[i/4 + j]);
            }

            let key = K::to_key(&chunk);
            
            let value = match index.get_mut(&key) {
                Some(v) => v,
                None => {
                    let v = Vec::new();
                    index.insert(key, v);
                    index.get_mut(&key).unwrap()
                }
            };

            value.push(i);
        }

        return IndexedPackedSequence { parent, index };
    }

    #[inline(always)]
    fn len(&self) -> usize {
        self.parent.len()
    }

    #[inline(always)]
    fn get(&self, idx: usize) -> char {
        self.parent.get(idx)
    }

    fn find_bounded(&self, pat: &impl PackedSequence, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        let pat_len = pat.len();

        if pat_len == 0 {
            panic!("Cannot search for an empty sequence, the return value is ambiguous");
        }

        // this check isn't really necessary because zero-length sequences will have empty indices,
        // and since it takes more than 0 cycles, it's getting axed.
        /*if self.len() == 0 {
            return None;
        }*/

        if C > pat_len {
            println!("Warning: chunk length is greater than the length of the item, falling back to slow implementation");
            return self.parent.find_bounded(pat, start, end);
        }

        let pat_packed = pat.get_packed();
        let key = K::to_key(pat_packed);

        if !self.index.contains_key(&key) {
            return None;
        }

        let start = match start {
            Some(s) => s,
            None => 0
        };
        let end = match end {
            Some(e) => e.min(self.len()),
            None => self.len()
        } - pat_len; // actual logic below is `if i + pat_len > end { break; }` but it is faster to
                     // just do the pat_len offset here, instead of in the inner loop

        let comparison_section_mask = ((0xffu16 << (8 - 2 * (pat_len % 4))) & 0xff) as u8;
        let pat_bytes_floor = pat_len / 4;

        // TODO: need variant of this method that returns them all in one go
        for &i in &self.index[&key] {
            if i < start {
                continue;
            }
            if i > end {
                break;
            }

            let variant = &self.parent.variants[i % 4];
            // Skip K::BYTES_NEEDED_FLOOR because they have already been checked against the index
            let comparison_section = &variant[i/4 + K::BYTES_NEEDED_FLOOR .. i/4 + pat_bytes_floor];//.to_owned();

            // check all the 'full' bytes
            if &pat_packed[K::BYTES_NEEDED_FLOOR..pat_bytes_floor] != comparison_section {
                continue; // mismatch, continue to next index
            }

            // if there's a final (partial) byte that has to be checked
            if comparison_section_mask != 0 {
                let comparison_byte = variant[i/4 + pat_bytes_floor] & comparison_section_mask;
                if comparison_byte != pat_packed[pat_bytes_floor] {
                    continue;
                }
            }

            return Some(i);
        }

        return None;
    }

    // NOTE: this should stay as much aligned with find_bounded as possible
    fn find_all(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        let pat_len = pat.len();

        if pat_len == 0 {
            panic!("Cannot search for an empty sequence, the return value is ambiguous");
        }

        if C > pat_len {
            println!("Warning: chunk length is greater than the length of the item, falling back to slow implementation");
            return self.parent.find_all(pat);
        }

        let pat_packed = pat.get_packed();
        let key = K::to_key(pat_packed);

        if !self.index.contains_key(&key) {
            return None;
        }

        let end = self.len() - pat_len;

        let comparison_section_mask = ((0xffu16 << (8 - 2 * (pat_len % 4))) & 0xff) as u8;
        let pat_bytes_floor = pat_len / 4;

        let mut out = Vec::new();

        for &i in &self.index[&key] {
            if i > end {
                break;
            }

            let variant = &self.parent.variants[i % 4];
            let comparison_section = &variant[i/4 + K::BYTES_NEEDED_FLOOR .. i/4 + pat_bytes_floor];

            // check all the 'full' bytes
            if &pat_packed[K::BYTES_NEEDED_FLOOR..pat_bytes_floor] != comparison_section {
                continue; // mismatch, continue to next index
            }

            // if there's a final (partial) byte that has to be checked
            if comparison_section_mask != 0 {
                let comparison_byte = variant[i/4 + pat_bytes_floor] & comparison_section_mask;
                if comparison_byte != pat_packed[pat_bytes_floor] {
                    continue;
                }
            }

            out.push(i);
        }

        return Some(out);
    }

    unsafe fn find_all_31mer(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        const PAT_LEN: usize = 31;

        // WARN: the caller must check this
        /*if C > PAT_LEN {
            println!("Warning: chunk length is greater than the length of the item, falling back to slow implementation");
            return self.parent.find_all(pat);
        }*/

        let pat_packed = pat.get_packed();
        let key = K::to_key(pat_packed);

        if !self.index.contains_key(&key) {
            return None;
        }

        let end = self.len() - PAT_LEN;

        const COMPARISON_SECTION_MASK: u8 = ((0xffu16 << (8 - 2 * (PAT_LEN % 4))) & 0xff) as u8;
        const PAT_BYTES_FLOOR: usize = PAT_LEN / 4;

        let mut out = Vec::new();

        for &i in &self.index[&key] {
            if i > end {
                break;
            }

            let variant = &self.parent.variants[i % 4];
            let comparison_section = &variant[i/4 + K::BYTES_NEEDED_FLOOR .. i/4 + PAT_BYTES_FLOOR];

            // check all the 'full' bytes
            if &pat_packed[K::BYTES_NEEDED_FLOOR..PAT_BYTES_FLOOR] != comparison_section {
                continue; // mismatch, continue to next index
            }

            // if there's a final (partial) byte that has to be checked
            if COMPARISON_SECTION_MASK != 0 {
                let comparison_byte = variant[i/4 + PAT_BYTES_FLOOR] & COMPARISON_SECTION_MASK;
                if comparison_byte != pat_packed[PAT_BYTES_FLOOR] {
                    continue;
                }
            }

            out.push(i);
        }

        return Some(out);
    }

    #[inline(always)]
    fn subsections(&self, chunk_length: usize) -> impl Iterator<Item = SimplePackedSequence> {
        self.parent.subsections(chunk_length)
    }

    #[inline(always)]
    fn get_packed<'a>(&'a self) -> &'a Vec<u8> {
        self.parent.get_packed()
    }
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub enum StandardIndexedPackedSequence {
    CS1(IndexedPackedSequence<u8, 1>),
    CS2(IndexedPackedSequence<u8, 2>),
    CS3(IndexedPackedSequence<u8, 3>),
    CS4(IndexedPackedSequence<u8, 4>),
    CS5(IndexedPackedSequence<u16, 5>),
    CS6(IndexedPackedSequence<u16, 6>),
    CS7(IndexedPackedSequence<u16, 7>),
    CS8(IndexedPackedSequence<u16, 8>),
    CS16(IndexedPackedSequence<u32, 16>),
    CS32(IndexedPackedSequence<u64, 32>)
}

impl StandardIndexedPackedSequence {
    pub fn new(seq: &str, chunk_length: usize) -> Option<StandardIndexedPackedSequence> {
        match chunk_length {
            1 => Some(StandardIndexedPackedSequence::CS1(indexed_packed_sequence!(seq, 1))),
            2 => Some(StandardIndexedPackedSequence::CS2(indexed_packed_sequence!(seq, 2))),
            3 => Some(StandardIndexedPackedSequence::CS3(indexed_packed_sequence!(seq, 3))),
            4 => Some(StandardIndexedPackedSequence::CS4(indexed_packed_sequence!(seq, 4))),
            5 => Some(StandardIndexedPackedSequence::CS5(indexed_packed_sequence!(seq, 5))),
            6 => Some(StandardIndexedPackedSequence::CS6(indexed_packed_sequence!(seq, 6))),
            7 => Some(StandardIndexedPackedSequence::CS7(indexed_packed_sequence!(seq, 7))),
            8 => Some(StandardIndexedPackedSequence::CS8(indexed_packed_sequence!(seq, 8))),
            16 => Some(StandardIndexedPackedSequence::CS16(indexed_packed_sequence!(seq, 16))),
            32 => Some(StandardIndexedPackedSequence::CS32(indexed_packed_sequence!(seq, 32))),
            _ => None
        }
    }

    pub fn chunk_length(&self) -> usize {
        match self {
            StandardIndexedPackedSequence::CS1(_) => 1,
            StandardIndexedPackedSequence::CS2(_) => 2,
            StandardIndexedPackedSequence::CS3(_) => 3,
            StandardIndexedPackedSequence::CS4(_) => 4,
            StandardIndexedPackedSequence::CS5(_) => 5,
            StandardIndexedPackedSequence::CS6(_) => 6,
            StandardIndexedPackedSequence::CS7(_) => 7,
            StandardIndexedPackedSequence::CS8(_) => 8,
            StandardIndexedPackedSequence::CS16(_) => 16,
            StandardIndexedPackedSequence::CS32(_) => 32,
        }
    }
}

impl PackedSequence for StandardIndexedPackedSequence {
    fn new(seq: &str) -> Self where Self: Sized {
        Self::CS8(IndexedPackedSequence::<u16, 8>::new(seq))
    }

    fn len(&self) -> usize {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.len()
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.len()
            },
        }
    }

    fn get(&self, idx: usize) -> char {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.get(idx)
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.get(idx)
            },
        }
    }

    fn find_bounded(&self, pat: &impl PackedSequence, start: Option<usize>, end: Option<usize>) -> Option<usize> {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.find_bounded(pat, start, end)
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.find_bounded(pat, start, end)
            },
        }
    }

    fn subsections(&self, chunk_length: usize) -> impl Iterator<Item = SimplePackedSequence> {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.subsections(chunk_length).collect::<Vec<_>>().into_iter()
            },
        }
    }

    fn get_packed<'a>(&'a self) -> &'a Vec<u8> {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.get_packed()
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.get_packed()
            },
        }
    }

    fn find_all(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.find_all(pat)
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.find_all(pat)
            },
        }
    }

    unsafe fn find_all_31mer(&self, pat: &impl PackedSequence) -> Option<Vec<usize>> {
        match self {
            StandardIndexedPackedSequence::CS1(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS2(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS3(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS4(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS5(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS6(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS7(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS8(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS16(v) => {
                v.find_all_31mer(pat)
            },
            StandardIndexedPackedSequence::CS32(v) => {
                v.find_all_31mer(pat)
            },
        }
    }
}

#[macro_export]
macro_rules! indexed_packed_sequence {
    ($seq:expr, 1) => {
        IndexedPackedSequence::<u8, 1>::new($seq)
    };
    ($seq:expr, 2) => {
        IndexedPackedSequence::<u8, 2>::new($seq)
    };
    ($seq:expr, 3) => {
        IndexedPackedSequence::<u8, 3>::new($seq)
    };
    ($seq:expr, 4) => {
        IndexedPackedSequence::<u8, 4>::new($seq)
    };
    ($seq:expr, 5) => {
        IndexedPackedSequence::<u16, 5>::new($seq)
    };
    ($seq:expr, 6) => {
        IndexedPackedSequence::<u16, 6>::new($seq)
    };
    ($seq:expr, 7) => {
        IndexedPackedSequence::<u16, 7>::new($seq)
    };
    ($seq:expr, 8) => {
        IndexedPackedSequence::<u16, 8>::new($seq)
    };
    ($seq:expr, 16) => {
        IndexedPackedSequence::<u32, 16>::new($seq)
    };
    ($seq:expr, 32) => {
        IndexedPackedSequence::<u64, 32>::new($seq)
    };
}

#[cfg(test)]
mod tests;
