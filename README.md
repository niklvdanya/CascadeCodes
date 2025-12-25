# CascadeCodes

Concatenated error-correcting codes in Rust.

## Features

- Hamming codes (7,4), (15,11), (31,26), etc.
- Reed-Solomon codes over GF(16) and GF(256)
- Block interleaver
- Two-level cascade codes
- Multi-level cascade codes (N levels)

## Build

```bash
cargo build --release
```

## Run demos

```bash
cargo run --bin cascade-cli demo all
cargo run --bin cascade-cli bench
```

## Usage

```rust
use cascade_codes::prelude::*;

// Simple cascade
let cascade = CascadeCode::new(
    ReedSolomonCode::new(255, 223, 8).unwrap(),
    HammingCode::hamming_7_4().unwrap(),
);
let encoded = cascade.encode(b"data").unwrap();
let decoded = cascade.decode(&encoded).unwrap();

// Multi-level cascade
let ml = MultiLevelCascade::new()
    .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
    .add_level_with_interleaver(
        ReedSolomonCode::new(15, 9, 4).unwrap(),
        Box::new(BlockInterleaver::new(16, 16)),
    )
    .add_level(HammingCode::hamming_7_4().unwrap());
let encoded = ml.encode(b"data").unwrap();
let decoded = ml.decode(&encoded).unwrap();
```
