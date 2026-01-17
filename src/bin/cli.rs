use std::env;
use cascade_codes::prelude::*;

fn main() {
    match env::args().nth(1).as_deref() {
        Some("demo") => run_demo(env::args().nth(2).as_deref().unwrap_or("all")),
        Some("bench") => run_bench(),
        Some("test") => run_tests(),
        _ => println!("Usage: cascade-cli [demo|bench|test] [hamming|rs|cascade|all]"),
    }
}

fn run_tests() {
    println!("\n=== Tests ===\n");
    
    let mut passed = 0;
    let mut failed = 0;
    
    print!("Test 1: MultiLevel 2-level (no errors)... ");
    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap());
    let msg = b"Hello, World!";
    let encoded = ml.encode(msg).unwrap();
    let decoded = ml.decode(&encoded).unwrap();
    if decoded == msg { println!("OK"); passed += 1; } 
    else { println!("FAIL"); failed += 1; }
    
    print!("Test 2: MultiLevel 3-level (no errors)... ");
    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 247, 8).unwrap());
    let msg = b"Triple cascade test!";
    let encoded = ml.encode(msg).unwrap();
    match ml.decode(&encoded) {
        Ok(d) if d == msg => { println!("OK"); passed += 1; }
        _ => { println!("FAIL"); failed += 1; }
    }
    
    print!("Test 3: MultiLevel with errors... ");
    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap());
    let msg = b"Error test";
    let encoded = ml.encode(msg).unwrap();
    let mut corrupted = encoded.clone();
    corrupted[10] ^= 0xFF;
    corrupted[20] ^= 0xFF;
    match ml.decode(&corrupted) {
        Ok(d) if &d[..] == msg => { println!("OK"); passed += 1; }
        _ => { println!("FAIL"); failed += 1; }
    }
    
    print!("Test 4: CascadeCode RS + RS... ");
    let cascade = CascadeCode::new(
        ReedSolomonCode::new(255, 223, 8).unwrap(),
        ReedSolomonCode::new(255, 239, 8).unwrap(),
    );
    let msg = b"Cascade RS+RS";
    let encoded = cascade.encode(msg).unwrap();
    let mut corrupted = encoded.clone();
    corrupted[50] ^= 0xFF;
    corrupted[100] ^= 0xFF;
    match cascade.decode(&corrupted) {
        Ok(d) if &d[..] == msg => { println!("OK"); passed += 1; }
        _ => { println!("FAIL"); failed += 1; }
    }
    
    print!("Test 5: CascadeCode RS + Hamming... ");
    let cascade = CascadeCode::with_converter(
        ReedSolomonCode::new(255, 223, 8).unwrap(),
        HammingCode::hamming_7_4().unwrap(),
        Box::new(BitConverter::for_hamming_7_4()),
    );
    let msg = b"RS+Hamming";
    let encoded = cascade.encode(msg).unwrap();
    let mut corrupted = encoded.clone();
    corrupted[10] ^= 0b00000001;
    corrupted[20] ^= 0b00000001;
    corrupted[30] ^= 0b00000001;
    match cascade.decode(&corrupted) {
        Ok(d) if &d[..] == msg => { println!("OK"); passed += 1; }
        _ => { println!("FAIL"); failed += 1; }
    }
    
    print!("Test 6: Long message (1KB)... ");
    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap());
    let msg: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();
    let encoded = ml.encode(&msg).unwrap();
    match ml.decode(&encoded) {
        Ok(d) if d == msg => { println!("OK"); passed += 1; }
        _ => { println!("FAIL"); failed += 1; }
    }
    
    println!("\nResults: {} passed, {} failed", passed, failed);
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
    println!("\n--- Hamming(7,4) ---");
    
    let code = HammingCode::hamming_7_4().unwrap();
    let data = vec![0b00001101];
    
    println!("Input: 0b{:08b}", data[0]);
    
    let encoded = code.encode_block(&data).unwrap();
    println!("Encoded: 0b{:08b}", encoded[0]);
    
    let mut corrupted = encoded.clone();
    corrupted[0] ^= 0b00000100;
    println!("Corrupted: 0b{:08b} (bit 2 flipped)", corrupted[0]);
    
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("Decoded: 0b{:08b}", decoded[0]);
    
    let ok = (decoded[0] & 0x0F) == (data[0] & 0x0F);
    println!("Result: {}", if ok { "CORRECTED" } else { "FAILED" });
}

fn demo_rs() {
    println!("\n--- Reed-Solomon(15,9) ---");
    
    let code = ReedSolomonCode::new(15, 9, 4).unwrap();
    let data: Vec<u8> = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
    
    println!("Input: {:?}", data);
    println!("Can correct: {} symbols", code.t());
    
    let encoded = code.encode_block(&data).unwrap();
    println!("Encoded: {:?}", encoded);
    
    let mut corrupted = encoded.clone();
    corrupted[0] ^= 0x0F;
    corrupted[5] ^= 0x0A;
    corrupted[10] ^= 0x05;
    println!("Corrupted: {:?}", corrupted);
    
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("Decoded: {:?}", decoded);
    println!("Result: {}", if decoded == data { "CORRECTED" } else { "FAILED" });

    println!("\n--- Reed-Solomon(255,223) ---");
    
    let code = ReedSolomonCode::new(255, 223, 8).unwrap();
    let data: Vec<u8> = (0..223).map(|i| i as u8).collect();
    
    println!("Input: [0, 1, 2, ... 222] (223 bytes)");
    println!("Can correct: {} symbols", code.t());
    
    let encoded = code.encode_block(&data).unwrap();
    let mut corrupted = encoded.clone();
    for i in 0..16 {
        corrupted[i * 10] ^= 0xFF;
    }
    println!("Corrupted[0..20]: {:?}", &corrupted[0..20]);
    
    let decoded = code.decode_block(&corrupted).unwrap();
    println!("Result: {}", if decoded == data { "CORRECTED" } else { "FAILED" });
}

fn demo_cascade() {
    println!("\n--- Cascade: RS(255,223) + RS(255,239) ---");

    let code = CascadeCode::new(
        ReedSolomonCode::new(255, 223, 8).unwrap(),
        ReedSolomonCode::new(255, 239, 8).unwrap(),
    );

    println!("Converter: {}", code.converter_info());
    println!("Rate: {:.3}", code.rate());

    let msg = b"Cascade!";
    println!("Original: {:?}", msg);

    let encoded = code.encode(msg).unwrap();
    println!("Encoded: {} bytes", encoded.len());
    println!("Encoded[0..32]: {:?}", &encoded[0..32.min(encoded.len())]);

    let mut corrupted = encoded.clone();
    corrupted[10] ^= 0xFF;
    corrupted[20] ^= 0xFF;
    corrupted[30] ^= 0xFF;
    println!("Corrupted[0..32]: {:?}", &corrupted[0..32.min(corrupted.len())]);

    match code.decode(&corrupted) {
        Ok(decoded) => {
            println!("Decoded: {:?}", decoded);
            println!("Result: {}", if &decoded[..] == msg { "CORRECTED" } else { "FAILED" });
        }
        Err(e) => println!("Error: {:?}", e),
    }

    println!("\n--- Cascade: RS(255,223) + Hamming(7,4) ---");

    let code = CascadeCode::with_converter(
        ReedSolomonCode::new(255, 223, 8).unwrap(),
        HammingCode::hamming_7_4().unwrap(),
        Box::new(BitConverter::for_hamming_7_4()),
    );

    println!("Converter: {}", code.converter_info());
    println!("Rate: {:.3}", code.rate());

    let msg = b"Test!";
    println!("Original: {:?}", msg);

    let encoded = code.encode(msg).unwrap();
    println!("Encoded: {} bytes", encoded.len());
    println!("Encoded[0..20]: {:?}", &encoded[0..20.min(encoded.len())]);

    let mut corrupted = encoded.clone();
    corrupted[10] ^= 0b00000001;
    corrupted[20] ^= 0b00000001;
    corrupted[30] ^= 0b00000001;
    println!("Corrupted[0..20]: {:?}", &corrupted[0..20.min(corrupted.len())]);
    println!("Errors: 3 single-bit errors at pos 10, 20, 30");

    match code.decode(&corrupted) {
        Ok(decoded) => {
            println!("Decoded: {:?}", decoded);
            println!("Result: {}", if &decoded[..] == msg { "CORRECTED" } else { "FAILED" });
        }
        Err(e) => println!("Error: {:?}", e),
    }

    println!("\n--- Multi-Level: RS + RS + RS ---");

    let ml = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 247, 8).unwrap());

    let msg = b"Triple!";
    println!("Original: {:?}", msg);
    println!("Levels: {}, Rate: {:.3}", ml.levels(), ml.rate());

    let encoded = ml.encode(msg).unwrap();
    println!("Encoded: {} bytes", encoded.len());

    let mut corrupted = encoded.clone();
    let positions = [100, 200, 300];
    for &pos in &positions {
        if pos < corrupted.len() { corrupted[pos] ^= 0xFF; }
    }
    println!("Error positions: {:?}", positions);
    println!("Corrupted[98..108]: {:?}", &corrupted[98..108.min(corrupted.len())]);

    match ml.decode(&corrupted) {
        Ok(decoded) => {
            println!("Decoded: {:?}", decoded);
            println!("Result: {}", if &decoded[..] == msg { "CORRECTED" } else { "FAILED" });
        }
        Err(e) => println!("Error: {:?}", e),
    }

    println!("\n--- Interleaver Demo: Burst Error ---");
    
    let rs = ReedSolomonCode::new(255, 223, 8).unwrap();
    println!("RS(255,223): t={} errors per block", rs.t());
    
    let msg: Vec<u8> = (0u8..=255).cycle().take(2000).collect();
    let num_blocks = (msg.len() + 222) / 223;
    println!("Message: {} bytes -> {} blocks", msg.len(), num_blocks);
    
    let encoded = rs.encode(&msg).unwrap();
    
    let burst = 20;
    println!("Burst: {} consecutive errors (> t=16)", burst);
    
    let mut corrupted = encoded.clone();
    for i in 0..burst {
        corrupted[100 + i] ^= 0xFF;
    }
    println!("Corrupted[98..125]: {:?}", &corrupted[98..125]);
    
    print!("Without interleaver: ");
    match rs.decode(&corrupted) {
        Ok(d) if d == msg => println!("OK"),
        _ => println!("FAIL"),
    }
    
    fn interleave(data: &[u8], depth: usize) -> Vec<u8> {
        let mut out = vec![0u8; data.len()];
        for (i, &b) in data.iter().enumerate() {
            let pos = (i % depth) * ((data.len() + depth - 1) / depth) + i / depth;
            if pos < out.len() { out[pos] = b; }
        }
        out
    }
    
    fn deinterleave(data: &[u8], depth: usize) -> Vec<u8> {
        let mut out = vec![0u8; data.len()];
        for i in 0..data.len() {
            let pos = (i % depth) * ((data.len() + depth - 1) / depth) + i / depth;
            if pos < data.len() { out[i] = data[pos]; }
        }
        out
    }
    
    let interleaved = interleave(&encoded, num_blocks);
    let mut corrupted_il = interleaved.clone();
    for i in 0..burst { corrupted_il[100 + i] ^= 0xFF; }
    let deinterleaved = deinterleave(&corrupted_il, num_blocks);
    
    print!("With interleaver:    ");
    match rs.decode(&deinterleaved) {
        Ok(d) if d == msg => println!("OK"),
        _ => println!("FAIL"),
    }
    
    println!("{} errors / {} blocks = {} per block", burst, num_blocks, burst / num_blocks);
}

fn run_bench() {
    println!("\n=== Benchmark (1KB, 100x) ===\n");
    
    let data: Vec<u8> = (0..1000).map(|i| (i % 256) as u8).collect();

    let code = HammingCode::hamming_7_4().unwrap();
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e);
    }
    println!("Hamming(7,4):     {:?}", start.elapsed());

    let code = ReedSolomonCode::new(255, 223, 8).unwrap();
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e);
    }
    println!("RS(255,223):      {:?}", start.elapsed());

    let code = CascadeCode::with_converter(
        ReedSolomonCode::new(255, 223, 8).unwrap(),
        HammingCode::hamming_7_4().unwrap(),
        Box::new(BitConverter::for_hamming_7_4()),
    );
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e);
    }
    println!("RS + Hamming:     {:?}", start.elapsed());

    let code = MultiLevelCascade::new()
        .add_level(ReedSolomonCode::new(255, 223, 8).unwrap())
        .add_level(ReedSolomonCode::new(255, 239, 8).unwrap());
    let start = std::time::Instant::now();
    for _ in 0..100 {
        let e = code.encode(&data).unwrap();
        let _ = code.decode(&e).unwrap();
    }
    println!("MultiLevel 2-lvl: {:?}", start.elapsed());
}
