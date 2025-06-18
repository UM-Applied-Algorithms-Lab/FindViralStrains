use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::collections::BTreeMap;
use csv::Writer;

#[derive(Debug)]
struct AlignmentStats {
    sample_name: String,
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
    let mut results = BTreeMap::new();

    // Process each subgraph directory
    for entry in fs::read_dir(input_dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_dir() {
            if let Some(dir_name) = path.file_name() {
                let dir_name = dir_name.to_string_lossy();
                if dir_name.starts_with("subgraph_") {
                    if let Some((sample_name, stats)) = process_subgraph_dir(&path)? {
                        results.insert(dir_name.to_string(), (sample_name, stats));
                    }
                }
            }
        }
    }

    // Write CSV output
    write_csv_output(output_path, &results)?;

    println!("Successfully processed {} subgraphs, output written to {}", 
        results.len(), 
        output_path.display());

    Ok(())
}

fn process_subgraph_dir(dir: &Path) -> std::io::Result<Option<(String, AlignmentStats)>> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name = file_name.to_string_lossy();
                if file_name.ends_with("_vs_ref.txt") && file_name.contains("1_of_1") {
                    // Extract sample name from filename (e.g., E1250_S84_L001 from E1250_S84_L001_1_of_1_vs_ref.txt)
                    let sample_name = file_name.split('_')
                        .take(3)
                        .collect::<Vec<_>>()
                        .join("_");
                    return Ok(Some((sample_name, parse_alignment_file(&path)?)));
                }
            }
        }
    }
    Ok(None)
}

fn parse_alignment_file(file_path: &Path) -> std::io::Result<AlignmentStats> {
    let file = fs::File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut stats = AlignmentStats {
        sample_name: String::new(),
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

fn write_csv_output(output_path: &Path, results: &BTreeMap<String, (String, AlignmentStats)>) -> std::io::Result<()> {
    let mut writer = Writer::from_path(output_path)?;

    // Write header
    writer.write_record(&[
        "Subgraph",
        "Sample",
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

    // Write data with one row per subgraph
    for (subgraph, (sample_name, stats)) in results {
        let alignment_length = stats.end_position - stats.start_position + 1;
        writer.write_record(&[
            subgraph,
            sample_name,
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
