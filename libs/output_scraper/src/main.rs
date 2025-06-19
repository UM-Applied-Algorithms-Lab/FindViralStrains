use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::collections::BTreeMap;
use csv::Writer;

#[derive(Debug)]
struct AlignmentStats {
    sample_name: String,
    subgraph_name: String,  // Added to track subgraph directory
    part_number: usize,
    total_parts: usize,
    length: usize,
    identity_pct: f64,
    identity_count: usize,
    gaps_pct: f64,
    gaps_count: usize,
    score: f64,
    start_position: usize,
    end_position: usize,
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_directory> <output_csv>", args[0]);
        std::process::exit(1);
    }

    let input_dir = Path::new(&args[1]);
    let output_path = Path::new(&args[2]);
    let mut results = Vec::new();

    // Process each sample directory
    for sample_entry in fs::read_dir(input_dir)? {
        let sample_entry = sample_entry?;
        let sample_path = sample_entry.path();
        
        if sample_path.is_dir() {
            let sample_name = sample_path.file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();

            // Process each subgraph directory in the sample directory
            for subgraph_entry in fs::read_dir(&sample_path)? {
                let subgraph_entry = subgraph_entry?;
                let subgraph_path = subgraph_entry.path();
                
                if subgraph_path.is_dir() {
                    if let Some(dir_name) = subgraph_path.file_name() {
                        let subgraph_name = dir_name.to_string_lossy().to_string();
                        if subgraph_name.starts_with("subgraph_") {
                            if let Some(stats_vec) = process_subgraph_dir(&subgraph_path, &sample_name, &subgraph_name)? {
                                results.extend(stats_vec);
                            }
                        }
                    }
                }
            }
            
            // Also check for files directly in the sample directory (like subgraph_0 might be missing)
            if let Some(stats_vec) = process_files_in_dir(&sample_path, &sample_name, "root")? {
                results.extend(stats_vec);
            }
        }
    }

    // Write CSV output
    write_csv_output(output_path, &results)?;

    println!("Successfully processed {} alignment files, output written to {}", 
        results.len(), 
        output_path.display());

    Ok(())
}

fn process_subgraph_dir(dir: &Path, sample_name: &str, subgraph_name: &str) -> std::io::Result<Option<Vec<AlignmentStats>>> {
    process_files_in_dir(dir, sample_name, subgraph_name)
}

fn process_files_in_dir(dir: &Path, sample_name: &str, subgraph_name: &str) -> std::io::Result<Option<Vec<AlignmentStats>>> {
    let mut stats_vec = Vec::new();
    
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name = file_name.to_string_lossy();
                if file_name.ends_with("_vs_ref.txt") {
                    // Extract the part numbers
                    let part_numbers = extract_part_numbers(&file_name);
                    
                    stats_vec.push(parse_alignment_file(
                        &path, 
                        sample_name.to_string(),
                        subgraph_name.to_string(),
                        part_numbers
                    )?);
                }
            }
        }
    }
    
    if stats_vec.is_empty() {
        Ok(None)
    } else {
        Ok(Some(stats_vec))
    }
}

fn extract_part_numbers(filename: &str) -> (usize, usize) {
    let parts: Vec<&str> = filename.split('_').collect();
    for i in 0..parts.len() {
        if parts[i] == "of" && i > 0 && i < parts.len() - 1 {
            if let (Ok(current), Ok(total)) = (
                parts[i-1].parse::<usize>(),
                parts[i+1].parse::<usize>(),
            ) {
                return (current, total);
            }
        }
    }
    (1, 1) // Default values if parsing fails
}

fn parse_alignment_file(
    file_path: &Path, 
    sample_name: String,
    subgraph_name: String,
    part_numbers: (usize, usize)
) -> std::io::Result<AlignmentStats> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut stats = AlignmentStats {
        sample_name,
        subgraph_name,
        part_number: part_numbers.0,
        total_parts: part_numbers.1,
        length: 0,
        identity_pct: 0.0,
        identity_count: 0,
        gaps_pct: 0.0,
        gaps_count: 0,
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
            stats.identity_count = parse_count(identity_str);
        } 
        else if line.starts_with("# Gaps: ") {
            let gaps_str = line[8..].trim();
            stats.gaps_pct = parse_percentage(gaps_str);
            stats.gaps_count = parse_count(gaps_str);
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

fn parse_count(s: &str) -> usize {
    s.split('/').next()
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0)
}

fn write_csv_output(output_path: &Path, results: &[AlignmentStats]) -> std::io::Result<()> {
    let mut writer = Writer::from_path(output_path)?;

    // Write header
    writer.write_record(&[
        "Sample",
        "Subgraph",
        "Part",
        "Total Parts",
        "Length",
        "Identity %",
        "Identity Count",
        "Gaps %",
        "Gaps Count",
        "Score",
        "Start Position",
        "End Position",
        "Alignment Length",
    ])?;

    // Write data
    for stats in results {
        let alignment_length = stats.end_position - stats.start_position + 1;
        writer.write_record(&[
            &stats.sample_name,
            &stats.subgraph_name,
            &stats.part_number.to_string(),
            &stats.total_parts.to_string(),
            &stats.length.to_string(),
            &format!("{:.1}", stats.identity_pct),
            &stats.identity_count.to_string(),
            &format!("{:.1}", stats.gaps_pct),
            &stats.gaps_count.to_string(),
            &format!("{:.1}", stats.score),
            &stats.start_position.to_string(),
            &stats.end_position.to_string(),
            &alignment_length.to_string(),
        ])?;
    }

    writer.flush()?;
    Ok(())
}
