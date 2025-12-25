use std::env;
use cascade_codes::prelude::*;

fn main() {
    match env::args().nth(1).as_deref() {
        Some("demo") => run_demo(env::args().nth(2).as_deref().unwrap_or("all")),
        Some("bench") => run_bench(),
        _ => println!("Usage: cascade-cli [demo|bench] [hamming|rs|cascade|all]"),
    }
}

fn run_demo(t: &str) {
    match t {
        "hamming" => demo_hamming(),
        "rs" => demo_rs(),
        "cascade" => demo_cascade(),
        "all" => {
            demo_hamming();
            demo_rs();
            demo_cascade();
        }
        _ => println!("Unknown demo: {}", t),
    }
}

fn demo_hamming() {
    println!("\n[Hamming(7,4)]");
    let code = HammingCode::hamming_7_4().unwrap();
    let data = vec![1, 0, 1, 1];
    let encoded = code.encode_block(&data).unwrap();
    let mut corrupted = encoded.clone();
    corrupted[2] ^= 1;
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("input:     {:?}", data);
    println!("encoded:   {:?}", encoded);
    println!("corrupted: {:?}", corrupted);
    println!("decoded:   {:?}", decoded);
    println!("status:    {}", if decoded == data { "OK" } else { "FAIL" });
}

fn demo_rs() {
    println!("\n[Reed-Solomon(15,9) over GF(16)]");
    let code = ReedSolomonCode::new(15, 9, 4).unwrap();
    let data: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    let encoded = code.encode_block(&data).unwrap();
    let mut corrupted = encoded.clone();
    corrupted[0] ^= 0x0F;
    corrupted[5] ^= 0x0A;
    corrupted[10] ^= 0x05;
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("input:   {:?}", data);
    println!("3 errors injected at positions 0, 5, 10");
    println!("decoded: {:?}", decoded);
    println!("status:  {}", if decoded == data { "OK" } else { "FAIL" });

    println!("\n[Reed-Solomon(255,223) over GF(256)]");
    let code = ReedSolomonCode::new(255, 223, 8).unwrap();
    let data: Vec<u8> = (0..223).map(|i| i as u8).collect();
    let encoded = code.encode_block(&data).unwrap();
    let mut corrupted = encoded.clone();
    for i in 0..16 {
        corrupted[i * 10] ^= 0xFF;
    }
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("input:   [0..222] (223 bytes)");
    println!("16 errors injected");
    println!("status:  {}", if decoded == data { "OK" } else { "FAIL" });
}

fn demo_cascade() {
    println!("\n[Cascade RS+Hamming]");
    let code = CascadeCode::new(
        ReedSolomonCode::new(15, 9, 4).unwrap(),
        HammingCode::hamming_7_4().unwrap(),
    );
    let msg = b"Hello!";
    let encoded = code.encode(msg).unwrap();
    let decoded = code.decode(&encoded).unwrap();
    println!("input:   {:?}", String::from_utf8_lossy(msg));
    println!("encoded: {} bytes", encoded.len());
    println!("decoded: {:?}", String::from_utf8_lossy(&decoded));
    println!("status:  {}", if decoded == msg { "OK" } else { "FAIL" });

    println!("\n[MultiLevel: RS(255,223) + RS(15,9) + Hamming]");
    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level_with_interleaver(
            ReedSolomonCode::new(15, 9, 4).unwrap(),
            Box::new(BlockInterleaver::new(16, 16)),
        )
        .add_level(HammingCode::hamming_7_4().unwrap());

    let msg = b"Multi-level cascade code test!";
    let encoded = ml.encode(msg).unwrap();
    let decoded = ml.decode(&encoded).unwrap();
    println!("levels:  {}", ml.levels());
    println!("rate:    {:.3}", ml.rate());
    println!("input:   {:?}", String::from_utf8_lossy(msg));
    println!("encoded: {} bytes", encoded.len());
    println!("decoded: {:?}", String::from_utf8_lossy(&decoded));
    println!("status:  {}", if decoded == msg { "OK" } else { "FAIL" });
}

fn run_bench() {
    println!("\n[Benchmark: 1KB, 100 iterations]");
    let data: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();

    let code = HammingCode::hamming_7_4().unwrap();
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e);
    }
    println!("Hamming(7,4): {:?}", start.elapsed());

    let code = ReedSolomonCode::new(15, 9, 4).unwrap();
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e);
    }
    println!("RS(15,9):     {:?}", start.elapsed());
}
