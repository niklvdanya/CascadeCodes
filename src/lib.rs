use std::{fmt, error::Error};
use std::convert::TryInto;

#[derive(Debug, Clone)]
pub enum CodeError {
    InvalidParams(String),
    EncodingErr(String),
    DecodingErr(String),
    TooManyErrors(usize),
    BadBlockSize(usize, usize),
}

impl fmt::Display for CodeError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result { write!(f, "{:?}", self) }
}

impl Error for CodeError {}

pub type CodeResult<T> = Result<T, CodeError>;

#[derive(Debug, Clone)]
pub struct GaloisField {
    pub size: usize,
    pub m: usize,
    exp_table: Vec<u32>,
    log_table: Vec<i32>,
}

impl GaloisField {
    pub fn new(m: usize) -> CodeResult<Self> {
        let poly = match m {
            4 => 0b10011,
            5 => 0b100101,
            8 => 0b100011101,
            _ => return Err(CodeError::InvalidParams(format!("m={} not supported", m))),
        };
        let size = 1 << m;
        let mut exp = vec![0u32; size * 2];
        let mut log = vec![-1i32; size];
        let mut x = 1u32;
        for i in 0..(size - 1) {
            exp[i] = x;
            exp[i + size - 1] = x;
            log[x as usize] = i as i32;
            x <<= 1;
            if x >= size as u32 {
                x ^= poly;
                x &= (size - 1) as u32;
            }
        }
        Ok(GaloisField { size, m, exp_table: exp, log_table: log })
    }

    pub fn add(&self, a: u32, b: u32) -> u32 { a ^ b }

    pub fn mul(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 || a as usize >= self.size || b as usize >= self.size {
            return 0;
        }
        let la = self.log_table[a as usize];
        let lb = self.log_table[b as usize];
        if la < 0 || lb < 0 { return 0; }
        self.exp_table[((la + lb) as usize) % (self.size - 1)]
    }

    pub fn div(&self, a: u32, b: u32) -> u32 {
        if a == 0 || b == 0 || a as usize >= self.size || b as usize >= self.size {
            return 0;
        }
        let la = self.log_table[a as usize];
        let lb = self.log_table[b as usize];
        if la < 0 || lb < 0 { return 0; }
        let idx = (la - lb + (self.size - 1) as i32) as usize % (self.size - 1);
        self.exp_table[idx]
    }

    pub fn exp(&self, i: usize) -> u32 { self.exp_table[i % (self.size - 1)] }

    pub fn poly_eval(&self, p: &[u32], x: u32) -> u32 {
        p.iter().rev().fold(0, |acc, &c| self.add(self.mul(acc, x), c))
    }
}

pub trait BlockCode: Send + Sync {
    fn n(&self) -> usize;
    fn k(&self) -> usize;
    fn d(&self) -> usize;
    fn t(&self) -> usize { (self.d() - 1) / 2 }
    fn rate(&self) -> f64 { self.k() as f64 / self.n() as f64 }
    fn encode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>>;
    fn decode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>>;

    fn encode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        let mut out = (data.len() as u32).to_be_bytes().to_vec();
        for chunk in data.chunks(self.k()) {
            let mut blk = chunk.to_vec();
            blk.resize(self.k(), 0);
            out.extend(self.encode_block(&blk)?);
        }
        Ok(out)
    }

    fn decode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        if data.len() < 4 {
            return Err(CodeError::DecodingErr("too short".into()));
        }
        let len = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;
        let mut out = Vec::new();
        for chunk in data[4..].chunks(self.n()) {
            let mut blk = chunk.to_vec();
            blk.resize(self.n(), 0);
            out.extend(self.decode_block(&blk)?);
        }
        out.truncate(len);
        Ok(out)
    }
}

#[derive(Debug, Clone)]
pub struct HammingCode {
    r: usize,
    nn: usize,
    kk: usize,
}

impl HammingCode {
    pub fn new(n: usize, k: usize) -> CodeResult<Self> {
        let r = n - k;
        if (1 << r) - 1 != n {
            return Err(CodeError::InvalidParams("bad params".into()));
        }
        Ok(HammingCode { r, nn: n, kk: k })
    }

    pub fn hamming_7_4() -> CodeResult<Self> { Self::new(7, 4) }
    pub fn hamming_15_11() -> CodeResult<Self> { Self::new(15, 11) }
    pub fn hamming_31_26() -> CodeResult<Self> { Self::new(31, 26) }
}

impl BlockCode for HammingCode {
    fn n(&self) -> usize { (self.nn + 7) / 8 }
    fn k(&self) -> usize { (self.kk + 7) / 8 }
    fn d(&self) -> usize { 3 }

    fn encode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        let expected_bytes = (self.kk + 7) / 8;
        if data.len() != expected_bytes {
            return Err(CodeError::BadBlockSize(data.len(), expected_bytes));
        }
        let mut input_bits = Vec::with_capacity(self.kk);
        for &byte in data {
            for bit_pos in 0..8 {
                if input_bits.len() < self.kk {
                    input_bits.push((byte >> bit_pos) & 1);
                }
            }
        }
        let mut cw_bits = vec![0u8; self.nn];
        let mut data_idx = 0;
        for i in 1..=self.nn {
            if i.count_ones() != 1 {
                cw_bits[i - 1] = input_bits[data_idx];
                data_idx += 1;
            }
        }
        for p in 0..self.r {
            let pos = 1 << p;
            let mut parity = 0u8;
            for i in 1..=self.nn {
                if i & pos != 0 {
                    parity ^= cw_bits[i - 1];
                }
            }
            cw_bits[pos - 1] = parity;
        }
        let output_bytes = (self.nn + 7) / 8;
        let mut output = vec![0u8; output_bytes];
        for (i, &bit) in cw_bits.iter().enumerate() {
            if bit != 0 {
                output[i / 8] |= 1 << (i % 8);
            }
        }
        Ok(output)
    }

    fn decode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        let expected_bytes = (self.nn + 7) / 8;
        if data.len() != expected_bytes {
            return Err(CodeError::BadBlockSize(data.len(), expected_bytes));
        }
        let mut cw_bits = Vec::with_capacity(self.nn);
        for &byte in data {
            for bit_pos in 0..8 {
                if cw_bits.len() < self.nn {
                    cw_bits.push((byte >> bit_pos) & 1);
                }
            }
        }
        let mut syndrome = 0usize;
        for p in 0..self.r {
            let pos = 1 << p;
            let mut parity = 0u8;
            for i in 1..=self.nn {
                if i & pos != 0 {
                    parity ^= cw_bits[i - 1];
                }
            }
            if parity != 0 {
                syndrome |= pos;
            }
        }
        if syndrome > 0 && syndrome <= self.nn {
            cw_bits[syndrome - 1] ^= 1;
        }
        let mut output_bits = Vec::new();
        for i in 1..=self.nn {
            if i.count_ones() != 1 {
                output_bits.push(cw_bits[i - 1]);
            }
        }
        let output_bytes = (self.kk + 7) / 8;
        let mut output = vec![0u8; output_bytes];
        for (i, &bit) in output_bits.iter().enumerate() {
            if bit != 0 {
                output[i / 8] |= 1 << (i % 8);
            }
        }
        Ok(output)
    }
}

#[derive(Debug, Clone)]
pub struct ReedSolomonCode {
    nn: usize,
    kk: usize,
    tt: usize,
    gf: GaloisField,
    gen: Vec<u32>,
}

impl ReedSolomonCode {
    pub fn new(n: usize, k: usize, m: usize) -> CodeResult<Self> {
        let gf = GaloisField::new(m)?;
        let tt = (n - k) / 2;
        let mut gen = vec![1u32];
        for i in 0..(2 * tt) {
            let root = gf.exp(i);
            let mut ng = vec![0u32; gen.len() + 1];
            for (j, &c) in gen.iter().enumerate() {
                ng[j + 1] = gf.add(ng[j + 1], c);
                ng[j] = gf.add(ng[j], gf.mul(c, root));
            }
            gen = ng;
        }
        Ok(ReedSolomonCode { nn: n, kk: k, tt, gf, gen })
    }

    fn syndromes(&self, r: &[u32]) -> Vec<u32> {
        (0..2 * self.tt).map(|i| self.gf.poly_eval(r, self.gf.exp(i))).collect()
    }

    fn berlekamp(&self, s: &[u32]) -> Vec<u32> {
        let n2t = s.len();
        let mut c = vec![0u32; n2t + 1];
        c[0] = 1;
        let mut b = vec![0u32; n2t + 1];
        b[0] = 1;
        let mut l = 0usize;
        let mut m = 1usize;
        let mut bb = 1u32;
        for n in 0..n2t {
            let mut d = s[n];
            for i in 1..=l {
                d = self.gf.add(d, self.gf.mul(c[i], s[n - i]));
            }
            if d == 0 {
                m += 1;
            } else {
                let t = c.clone();
                let coef = self.gf.div(d, bb);
                for i in 0..n2t {
                    if i + m < c.len() && b[i] != 0 {
                        c[i + m] = self.gf.add(c[i + m], self.gf.mul(coef, b[i]));
                    }
                }
                if 2 * l <= n {
                    l = n + 1 - l;
                    b = t;
                    bb = d;
                    m = 1;
                } else {
                    m += 1;
                }
            }
        }
        c.truncate(l + 1);
        c
    }

    fn chien(&self, sigma: &[u32]) -> Vec<usize> {
        let mut pos = Vec::new();
        for i in 0..self.nn {
            let x = if i == 0 { 1 } else {
                self.gf.exp((self.gf.size - 1 - i) % (self.gf.size - 1))
            };
            let mut sum = 0u32;
            let mut xpow = 1u32;
            for &c in sigma.iter() {
                sum = self.gf.add(sum, self.gf.mul(c, xpow));
                xpow = self.gf.mul(xpow, x);
            }
            if sum == 0 {
                pos.push(i);
            }
        }
        pos
    }

    fn forney(&self, s: &[u32], sigma: &[u32], pos: &[usize]) -> Vec<u32> {
        let n2t = 2 * self.tt;
        let mut omega = vec![0u32; n2t];
        for i in 0..n2t {
            for j in 0..sigma.len().min(i + 1) {
                omega[i] = self.gf.add(omega[i], self.gf.mul(sigma[j], s[i - j]));
            }
        }
        let mut sigma_d = vec![0u32; sigma.len()];
        for i in (1..sigma.len()).step_by(2) {
            sigma_d[i - 1] = sigma[i];
        }
        pos.iter().map(|&p| {
            let x_p = self.gf.exp(p);
            let x_p_inv = if p == 0 { 1 } else {
                self.gf.exp((self.gf.size - 1 - p % (self.gf.size - 1)) % (self.gf.size - 1))
            };
            let mut ov = 0u32;
            let mut xpow = 1u32;
            for &c in omega.iter() {
                ov = self.gf.add(ov, self.gf.mul(c, xpow));
                xpow = self.gf.mul(xpow, x_p_inv);
            }
            let mut sv = 0u32;
            xpow = 1u32;
            for &c in sigma_d.iter() {
                sv = self.gf.add(sv, self.gf.mul(c, xpow));
                xpow = self.gf.mul(xpow, x_p_inv);
            }
            if sv != 0 { self.gf.mul(x_p, self.gf.div(ov, sv)) } else { 0 }
        }).collect()
    }
}

impl BlockCode for ReedSolomonCode {
    fn n(&self) -> usize { self.nn }
    fn k(&self) -> usize { self.kk }
    fn d(&self) -> usize { 2 * self.tt + 1 }

    fn encode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        if data.len() != self.kk {
            return Err(CodeError::BadBlockSize(data.len(), self.kk));
        }
        let red = self.nn - self.kk;
        let mut msg: Vec<u32> = data.iter().map(|&b| b as u32).collect();
        msg.resize(self.nn, 0);
        msg.rotate_right(red);
        let mut tmp = msg.clone();
        for i in (red..self.nn).rev() {
            if tmp[i] != 0 {
                let c = tmp[i];
                for (j, &g) in self.gen.iter().enumerate() {
                    tmp[i - (self.gen.len() - 1 - j)] =
                        self.gf.add(tmp[i - (self.gen.len() - 1 - j)], self.gf.mul(c, g));
                }
            }
        }
        let mut cw = vec![0u8; self.nn];
        for i in 0..red {
            cw[i] = tmp[i] as u8;
        }
        for i in 0..self.kk {
            cw[red + i] = data[i];
        }
        Ok(cw)
    }

    fn decode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        if data.len() != self.nn {
            return Err(CodeError::BadBlockSize(data.len(), self.nn));
        }
        let r: Vec<u32> = data.iter().map(|&b| b as u32).collect();
        let s = self.syndromes(&r);
        if s.iter().all(|&x| x == 0) {
            return Ok(data[self.nn - self.kk..].to_vec());
        }
        let sigma = self.berlekamp(&s);
        let pos = self.chien(&sigma);
        if pos.len() != sigma.len() - 1 {
            return Err(CodeError::TooManyErrors(pos.len()));
        }
        let mag = self.forney(&s, &sigma, &pos);
        let mut cor = r;
        for (i, &p) in pos.iter().enumerate() {
            cor[p] = self.gf.add(cor[p], mag[i]);
        }
        Ok(cor[self.nn - self.kk..].iter().map(|&x| x as u8).collect())
    }
}

pub trait Interleaver: Send + Sync {
    fn interleave(&self, data: &[u8]) -> Vec<u8>;
    fn deinterleave(&self, data: &[u8]) -> Vec<u8>;
}

#[derive(Debug, Clone)]
pub struct BlockInterleaver {
    rows: usize,
    cols: usize,
}

impl BlockInterleaver {
    pub fn new(rows: usize, cols: usize) -> Self {
        BlockInterleaver { rows, cols }
    }
}

impl Interleaver for BlockInterleaver {
    fn interleave(&self, data: &[u8]) -> Vec<u8> {
        let bs = self.rows * self.cols;
        data.chunks(bs).flat_map(|chunk| {
            let mut p = chunk.to_vec();
            p.resize(bs, 0);
            let mut out = vec![0u8; bs];
            for i in 0..self.rows {
                for j in 0..self.cols {
                    out[j * self.rows + i] = p[i * self.cols + j];
                }
            }
            out[..chunk.len().min(bs)].to_vec()
        }).collect()
    }

    fn deinterleave(&self, data: &[u8]) -> Vec<u8> {
        let bs = self.rows * self.cols;
        data.chunks(bs).flat_map(|chunk| {
            let mut p = chunk.to_vec();
            p.resize(bs, 0);
            let mut out = vec![0u8; bs];
            for i in 0..self.rows {
                for j in 0..self.cols {
                    out[i * self.cols + j] = p[j * self.rows + i];
                }
            }
            out[..chunk.len().min(bs)].to_vec()
        }).collect()
    }
}

pub trait LevelConverter: Send + Sync {
    fn convert_down(&self, data: &[u8]) -> Vec<u8>;
    fn convert_up(&self, data: &[u8]) -> Vec<u8>;
    fn info(&self) -> String { "converter".to_string() }
}

#[derive(Debug, Clone)]
pub struct IdentityConverter;

impl LevelConverter for IdentityConverter {
    fn convert_down(&self, data: &[u8]) -> Vec<u8> { data.to_vec() }
    fn convert_up(&self, data: &[u8]) -> Vec<u8> { data.to_vec() }
    fn info(&self) -> String { "identity".to_string() }
}

#[derive(Debug, Clone)]
pub struct BitConverter {
    bits_per_block: usize,
}

impl BitConverter {
    pub fn new(bits_per_block: usize) -> Self {
        BitConverter { bits_per_block }
    }
    
    pub fn for_hamming_7_4() -> Self { BitConverter { bits_per_block: 4 } }
    pub fn for_hamming_15_11() -> Self { BitConverter { bits_per_block: 11 } }
}

impl LevelConverter for BitConverter {
    fn convert_down(&self, data: &[u8]) -> Vec<u8> {
        let mut bits = Vec::new();
        for &byte in data {
            for i in 0..8 {
                bits.push((byte >> i) & 1);
            }
        }
        while bits.len() % self.bits_per_block != 0 {
            bits.push(0);
        }
        let bytes_per_block = (self.bits_per_block + 7) / 8;
        let mut result = Vec::new();
        for chunk in bits.chunks(self.bits_per_block) {
            let mut block = vec![0u8; bytes_per_block];
            for (i, &bit) in chunk.iter().enumerate() {
                if bit != 0 {
                    block[i / 8] |= 1 << (i % 8);
                }
            }
            result.extend(block);
        }
        result
    }
    
    fn convert_up(&self, data: &[u8]) -> Vec<u8> {
        let bytes_per_block = (self.bits_per_block + 7) / 8;
        let mut bits = Vec::new();
        for chunk in data.chunks(bytes_per_block) {
            let mut count = 0;
            for &byte in chunk {
                for i in 0..8 {
                    if count < self.bits_per_block {
                        bits.push((byte >> i) & 1);
                        count += 1;
                    }
                }
            }
        }
        let mut result = Vec::new();
        for chunk in bits.chunks(8) {
            let mut byte = 0u8;
            for (i, &bit) in chunk.iter().enumerate() {
                if bit != 0 {
                    byte |= 1 << i;
                }
            }
            result.push(byte);
        }
        result
    }
    
    fn info(&self) -> String {
        format!("bit converter ({} bits/block)", self.bits_per_block)
    }
}

pub struct CascadeCode<O: BlockCode, I: BlockCode> {
    pub outer: O,
    pub inner: I,
    interleaver: Option<Box<dyn Interleaver>>,
    converter: Box<dyn LevelConverter>,
}

impl<O: BlockCode, I: BlockCode> CascadeCode<O, I> {
    pub fn new(outer: O, inner: I) -> Self {
        CascadeCode { outer, inner, interleaver: None, converter: Box::new(IdentityConverter) }
    }
    
    pub fn with_interleaver(outer: O, inner: I, il: Box<dyn Interleaver>) -> Self {
        CascadeCode { outer, inner, interleaver: Some(il), converter: Box::new(IdentityConverter) }
    }
    
    pub fn with_converter(outer: O, inner: I, conv: Box<dyn LevelConverter>) -> Self {
        CascadeCode { outer, inner, interleaver: None, converter: conv }
    }
    
    pub fn with_all(outer: O, inner: I, il: Box<dyn Interleaver>, conv: Box<dyn LevelConverter>) -> Self {
        CascadeCode { outer, inner, interleaver: Some(il), converter: conv }
    }
    
    pub fn rate(&self) -> f64 { self.outer.rate() * self.inner.rate() }
    
    pub fn converter_info(&self) -> String { self.converter.info() }

    pub fn encode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        let outer_enc = self.outer.encode(data)?;
        let interleaved = self.interleaver.as_ref()
            .map_or(outer_enc.clone(), |il| il.interleave(&outer_enc));
        let converted = self.converter.convert_down(&interleaved);
        let inner_enc = self.inner.encode(&converted)?;
        Ok(inner_enc)
    }

    pub fn decode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        let inner_dec = self.inner.decode(data)?;
        let converted = self.converter.convert_up(&inner_dec);
        let deinterleaved = self.interleaver.as_ref()
            .map_or(converted.clone(), |il| il.deinterleave(&converted));
        let outer_dec = self.outer.decode(&deinterleaved)?;
        Ok(outer_dec)
    }
}

pub struct DynBlockCode(Box<dyn BlockCode>);

impl DynBlockCode {
    pub fn new<C: BlockCode + 'static>(code: C) -> Self { DynBlockCode(Box::new(code)) }
}

impl BlockCode for DynBlockCode {
    fn n(&self) -> usize { self.0.n() }
    fn k(&self) -> usize { self.0.k() }
    fn d(&self) -> usize { self.0.d() }
    fn encode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> { self.0.encode_block(data) }
    fn decode_block(&self, data: &[u8]) -> CodeResult<Vec<u8>> { self.0.decode_block(data) }
    fn encode(&self, data: &[u8]) -> CodeResult<Vec<u8>> { self.0.encode(data) }
    fn decode(&self, data: &[u8]) -> CodeResult<Vec<u8>> { self.0.decode(data) }
}

pub struct MultiLevelCascade {
    codes: Vec<DynBlockCode>,
    interleavers: Vec<Option<Box<dyn Interleaver>>>,
}

impl MultiLevelCascade {
    pub fn new() -> Self {
        MultiLevelCascade { codes: Vec::new(), interleavers: Vec::new() }
    }

    pub fn add_level<C: BlockCode + 'static>(mut self, code: C) -> Self {
        self.codes.push(DynBlockCode::new(code));
        self.interleavers.push(None);
        self
    }

    pub fn add_level_with_interleaver<C: BlockCode + 'static>(mut self, code: C, il: Box<dyn Interleaver>) -> Self {
        self.codes.push(DynBlockCode::new(code));
        self.interleavers.push(Some(il));
        self
    }

    pub fn levels(&self) -> usize { self.codes.len() }

    pub fn rate(&self) -> f64 { self.codes.iter().map(|c| c.rate()).product() }

    pub fn encode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        if self.codes.is_empty() {
            return Ok(data.to_vec());
        }
        let mut current = data.to_vec();
        let mut lengths: Vec<u32> = Vec::new();
        for (i, code) in self.codes.iter().enumerate() {
            lengths.push(current.len() as u32);
            let mut encoded = Vec::new();
            for chunk in current.chunks(code.k()) {
                let mut blk = chunk.to_vec();
                blk.resize(code.k(), 0);
                encoded.extend(code.encode_block(&blk)?);
            }
            if let Some(il) = &self.interleavers[i] {
                encoded = il.interleave(&encoded);
            }
            current = encoded;
        }
        let mut result = Vec::new();
        result.extend((self.codes.len() as u32).to_be_bytes());
        for len in &lengths {
            result.extend(len.to_be_bytes());
        }
        result.extend(current);
        Ok(result)
    }

    pub fn decode(&self, data: &[u8]) -> CodeResult<Vec<u8>> {
        if self.codes.is_empty() {
            return Ok(data.to_vec());
        }
        let header_size = 4 + 4 * self.codes.len();
        if data.len() < header_size {
            return Err(CodeError::DecodingErr(format!("data too short: {} < {}", data.len(), header_size)));
        }
        let num_levels = u32::from_be_bytes(data[0..4].try_into().unwrap()) as usize;
        if num_levels != self.codes.len() {
            return Err(CodeError::DecodingErr(format!("level mismatch: expected {}, got {}", self.codes.len(), num_levels)));
        }
        let mut lengths = Vec::new();
        for i in 0..num_levels {
            let offset = 4 + i * 4;
            let len = u32::from_be_bytes(data[offset..offset+4].try_into().unwrap()) as usize;
            lengths.push(len);
        }
        let mut current = data[header_size..].to_vec();
        for (i, code) in self.codes.iter().enumerate().rev() {
            if let Some(il) = &self.interleavers[i] {
                current = il.deinterleave(&current);
            }
            let mut decoded = Vec::new();
            for chunk in current.chunks(code.n()) {
                let mut blk = chunk.to_vec();
                blk.resize(code.n(), 0);
                decoded.extend(code.decode_block(&blk)?);
            }
            decoded.truncate(lengths[i]);
            current = decoded;
        }
        Ok(current)
    }
}

impl Default for MultiLevelCascade {
    fn default() -> Self { Self::new() }
}

pub mod prelude {
    pub use super::{
        BlockCode, BlockInterleaver, CascadeCode, CodeError, CodeResult,
        DynBlockCode, GaloisField, HammingCode, Interleaver, MultiLevelCascade,
        ReedSolomonCode, LevelConverter, IdentityConverter, BitConverter,
    };
}
