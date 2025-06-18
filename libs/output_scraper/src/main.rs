use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::collections::BTreeMap; // Changed from HashMap for ordered output

#[derive(Debug)]
struct AlignmentStats {
    file_name: String,
    length: usize,
    identity_pct: f64,
    gaps_pct: f64,
    score: f64,
    start_position: usize,
    end_position: usize,
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_directory> <output_file>", args[0]);
        std::process::exit(1);
    }

    let input_dir = Path::new(&args[1]);
    let output_path = Path::new(&args[2]);
    let mut results = BTreeMap::new(); // Key: subgraph name, Value: Vec of stats

    // Walk through the directory structure
    for entry in fs::read_dir(input_dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_dir() {
            if let Some(dir_name) = path.file_name() {
                let dir_name = dir_name.to_string_lossy();
                if dir_name.starts_with("subgraph_") {
                    process_subgraph_dir(&path, &dir_name, &mut results)?;
                }
            }
        }
    }

    // Write formatted results
    write_formatted_output(output_path, &results)?;

    println!("Successfully processed {} subgraphs, output written to {}", 
        results.len(), 
        output_path.display());

    Ok(())
}

fn process_subgraph_dir(dir: &Path, subgraph: &str, results: &mut BTreeMap<String, Vec<AlignmentStats>>) -> std::io::Result<()> {
    let mut subgraph_results = Vec::new();
    
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name = file_name.to_string_lossy();
                if file_name.ends_with("_vs_ref.txt") && file_name.contains("1_of_1") {
                    if let Ok(stats) = parse_alignment_file(&path) {
                        subgraph_results.push(stats);
                    }
                }
            }
        }
    }
    
    if !subgraph_results.is_empty() {
        results.insert(subgraph.to_string(), subgraph_results);
    }
    Ok(())
}

fn parse_alignment_file(file_path: &Path) -> std::io::Result<AlignmentStats> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);
    
    let file_name = file_path.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string();

    let mut stats = AlignmentStats {
        file_name,
        length: 0,
        identity_pct: 0.0,
        gaps_pct: 0.0,
        score: 0.0,
        start_position: 0,
        end_position: 0,
    };

    for line in reader.lines() {
        let line = line?;
        
        if line.starts_with("# Length: ") {
            stats.length = line[10..].trim().parse().unwrap_or(0);
        } 
        else if line.starts_with("# Identity: ") {
            let identity_str = line[12..].trim();
            stats.identity_pct = parse_percentage(identity_str);
        } 
        else if line.starts_with("# Gaps: ") {
            let gaps_str = line[8..].trim();
            stats.gaps_pct = parse_percentage(gaps_str);
        } 
        else if line.starts_with("# Score: ") {
            stats.score = line[9..].trim().parse().unwrap_or(0.0);
        } 
        else if line.starts_with("NC_045512.2") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                if stats.start_position == 0 {
                    stats.start_position = parts[1].parse().unwrap_or(0);
                }
                stats.end_position = parts.last().and_then(|s| s.parse().ok()).unwrap_or(0);
            }
        }
    }

    Ok(stats)
}

fn parse_percentage(s: &str) -> f64 {
    s.split('(').nth(1)
        .and_then(|s| s.split('%').next())
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0)
}

fn write_formatted_output(output_path: &Path, results: &BTreeMap<String, Vec<AlignmentStats>>) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;
    
    for (subgraph, stats_vec) in results {
        writeln!(file, "╔══════════════════════════════════════╗")?;
        writeln!(file, "║ Subgraph: {:<26} ║", subgraph)?;
        writeln!(file, "╠══════════════════════════════════════╣")?;
        
        for stats in stats_vec {
            writeln!(file, "║ File: {:<30} ║", stats.file_name)?;
            writeln!(file, "║   Length: {:<26} ║", stats.length)?;
            writeln!(file, "║   Identity: {:>5.1}% {:<18} ║", 
                stats.identity_pct, 
                format!("({}/{})", (stats.identity_pct/100.0 * stats.length as f64) as usize, stats.length))?;
            writeln!(file, "║   Gaps: {:>5.1}% {:<20} ║", 
                stats.gaps_pct,
                format!("({}/{})", (stats.gaps_pct/100.0 * stats.length as f64) as usize, stats.length))?;
            writeln!(file, "║   Score: {:<26.1} ║", stats.score)?;
            writeln!(file, "║   Positions: {}-{:<18} ║", stats.start_position, stats.end_position)?;
            writeln!(file, "╠──────────────────────────────────────╣")?;
        }
    }
    
    writeln!(file, "╚══════════════════════════════════════╝")?;
    Ok(())
}
