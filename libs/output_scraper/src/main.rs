use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::collections::{HashMap, HashSet};
use csv::Writer;

#[derive(Debug)]
struct AlignmentStats {
    sample_name: String,
    subgraph_name: String,
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
    runtime: f64,
    objective_value: f64,
    nodes: usize,
    edges: usize,
    sources: usize,
    sinks: usize,
    total_flow: f64, 
    explained_flow: f64, 
    weight: f64
}

#[derive(Debug, Clone)]
struct DecompStats {
    runtime: f64,
    objective_value: f64,
    total_flow: f64, 
    explained_flow: f64, 
    weight: f64
}

#[derive(Debug, Default)]
struct GraphData {
    nodes: HashSet<usize>,
    edges: usize,
    sources: usize,
    sinks: usize,
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input_directory> <output_csv>", args[0]);
        std::process::exit(1);
    }

    let input_dir = Path::new(&args[1]);
    let decomp_dir = input_dir.join("../decomp_results");
    let graphs_dir = input_dir.join("../graphs");
    let output_path = Path::new(&args[2]);
    let mut results = Vec::new();

    println!("Starting processing with:");
    println!("- Input directory: {}", input_dir.display());
    println!("- Decomp directory: {}", decomp_dir.display());
    println!("- Graphs directory: {}", graphs_dir.display());
    println!("- Output CSV: {}", output_path.display());

    // First collect all decomp stats in a lookup table
    println!("\nBuilding decomp stats map...");
    let decomp_stats_map = build_decomp_stats_map(&decomp_dir)?;
    println!("Found {} decomp results", decomp_stats_map.len());

    // Process each sample directory
    println!("\nProcessing sample directories...");
    for sample_entry in fs::read_dir(input_dir)? {
        let sample_entry = sample_entry?;
        let sample_path = sample_entry.path();
        
        if sample_path.is_dir() {
            let sample_name = sample_path.file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();

            println!("\nProcessing sample: {}", sample_name);

            // Process each subgraph directory in the sample directory
            println!("Processing subgraph directories...");
            for subgraph_entry in fs::read_dir(&sample_path)? {
                let subgraph_entry = subgraph_entry?;
                let subgraph_path = subgraph_entry.path();
                
                if subgraph_path.is_dir() {
                    if let Some(dir_name) = subgraph_path.file_name() {
                        let subgraph_name = dir_name.to_string_lossy().to_string();
                        if subgraph_name.starts_with("subgraph_") {
                            println!("  Processing subgraph: {}", subgraph_name);
                            if let Some(mut stats_vec) = process_subgraph_dir(&subgraph_path, &sample_name, &subgraph_name, &graphs_dir)? {
                                println!("    Found {} alignment files", stats_vec.len());
                                add_decomp_stats(&decomp_stats_map, &mut stats_vec);
                                results.extend(stats_vec);
                            }
                        }
                    }
                }
            }
            
            // Also check for files directly in the sample directory
            println!("Checking for root-level alignment files...");
            if let Some(mut stats_vec) = process_files_in_dir(&sample_path, &sample_name, "root", &graphs_dir)? {
                println!("  Found {} root-level alignment files", stats_vec.len());
                add_decomp_stats(&decomp_stats_map, &mut stats_vec);
                results.extend(stats_vec);
            }
        }
    }

    // Sort results by sample, then subgraph, then total parts, then part number
    results.sort_by(|a, b| {
        a.sample_name.cmp(&b.sample_name)
            .then(a.subgraph_name.cmp(&b.subgraph_name))
            .then(a.total_parts.cmp(&b.total_parts))
            .then(a.part_number.cmp(&b.part_number))
    });

    // Write CSV output
    println!("\nWriting output to {}...", output_path.display());
    write_csv_output(output_path, &results)?;

    println!("\nSuccessfully processed {} alignment files, output written to {}", 
        results.len(), 
        output_path.display());

    Ok(())
}



fn parse_graph_file(file_path: &Path) -> std::io::Result<GraphData> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut graph_data = GraphData::default();

    // Skip header line
    let mut lines = reader.lines().skip(1);

    while let Some(Ok(line)) = lines.next() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let (Ok(from_node), Ok(to_node)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                graph_data.nodes.insert(from_node);
                graph_data.nodes.insert(to_node);
                graph_data.edges += 1;
                
                // Count sources (edges from node 0)
                if from_node == 0 {
                    graph_data.sources += 1;
                }
                // Count sinks (edges to node 1)
                if to_node == 1 {
                    graph_data.sinks += 1;
                }
            }
        }
    }

    Ok(graph_data)
}

fn process_subgraph_dir(dir: &Path, sample_name: &str, subgraph_name: &str, graphs_dir: &Path) -> std::io::Result<Option<Vec<AlignmentStats>>> {
    process_files_in_dir(dir, sample_name, subgraph_name, graphs_dir)
}

fn process_files_in_dir(dir: &Path, sample_name: &str, subgraph_name: &str, graphs_dir: &Path) -> std::io::Result<Option<Vec<AlignmentStats>>> {
    let mut stats_vec = Vec::new();
    
    println!("    Scanning directory: {}", dir.display());
    
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name = file_name.to_string_lossy();
                if file_name.ends_with("_vs_ref.txt") {
                    println!("      Found alignment file: {}", file_name);
                    let part_numbers = extract_part_numbers(&file_name);
                    
                    let mut stats = parse_alignment_file(
                        &path, 
                        sample_name.to_string(),
                        subgraph_name.to_string(),
                        part_numbers
                    )?;
                    
                    // Find and parse graph file in the new format: <sample>.super_<num>.dbg
                    let subgraph_num = subgraph_name.trim_start_matches("subgraph_");
                    let graph_file_path = graphs_dir.join(format!("{}.super_{}.dbg", sample_name, subgraph_num));
                    
                    if graph_file_path.exists() {
                        println!("        Parsing graph file: {}", graph_file_path.display());
                        let graph_data = parse_graph_file(&graph_file_path)?;
                        stats.nodes = graph_data.nodes.len();
                        stats.edges = graph_data.edges;
                        stats.sources = graph_data.sources;
                        stats.sinks = graph_data.sinks;
                        
                        println!("          Nodes: {}, Edges: {}, Sources (from 0): {}, Sinks (to 1): {}", 
                            stats.nodes, stats.edges, stats.sources, stats.sinks);
                    } else {
                        println!("        Graph file not found: {}", graph_file_path.display());
                    }
                    
                    println!("        Part {}/{}: length={}, identity={:.1}%, gaps={:.1}%", 
                        stats.part_number, stats.total_parts, stats.length, 
                        stats.identity_pct, stats.gaps_pct);
                    
                    stats_vec.push(stats);
                }
            }
        }
    }
    
    if stats_vec.is_empty() {
        println!("      No alignment files found");
        Ok(None)
    } else {
        Ok(Some(stats_vec))
    }
}

fn build_decomp_stats_map(decomp_dir: &Path) -> std::io::Result<HashMap<(String, String, usize), DecompStats>> {
    let mut map = HashMap::new();
    
    println!("Scanning decomp directory: {}", decomp_dir.display());
    
    for entry in fs::read_dir(decomp_dir)? {
        let entry = entry?;
        let path = entry.path();
        
        if let Some(file_name) = path.file_name() {
            let file_name = file_name.to_string_lossy();
            if file_name.ends_with(".paths") {
                println!("  Found decomp file: {}", file_name);
                if let Some((sample_name, subgraph_name, total_parts)) = parse_decomp_filename(&file_name) {
                    println!("    Sample: {}, Subgraph: {}, Total Parts: {}", 
                        sample_name, subgraph_name, total_parts);

                    if let Some(stats) = parse_decomp_file(&path)? {
                        println!("    Runtime: {:.4}s, Objective: {:.6}, Total Flow: {:.6}, Explained Flow: {:.6}, Weight: {:.6}", 
                            stats.runtime, stats.objective_value, stats.total_flow, stats.explained_flow, stats.weight);
                        map.insert((sample_name, subgraph_name, total_parts), stats);
                    }
                }
            }
        }
    }
    Ok(map)
}

fn parse_decomp_filename(filename: &str) -> Option<(String, String, usize)> {
    let parts: Vec<&str> = filename.split('_').collect();
    if parts.len() >= 4 {
        let sample_end = parts.len() - 3;
        let sample_name = parts[..sample_end].join("_");
        let subgraph_name = format!("{}_{}", parts[sample_end], parts[sample_end + 1]);
        
        // Extract total parts from filename (assuming format like "XXX_YYY_Z.paths")
        let total_parts = parts.last()
            .and_then(|s| s.split('.').next())
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);
            
        return Some((sample_name, subgraph_name, total_parts));
    }
    None
}

fn add_decomp_stats(
    decomp_stats_map: &HashMap<(String, String, usize), DecompStats>,
    stats_vec: &mut Vec<AlignmentStats>
) {
    for stat in stats_vec {
        let key = (
            stat.sample_name.clone(), 
            stat.subgraph_name.clone(),
            stat.total_parts
        );
        if let Some(decomp_stats) = decomp_stats_map.get(&key) {
            stat.runtime = decomp_stats.runtime;
            stat.objective_value = decomp_stats.objective_value;
            stat.total_flow = decomp_stats.total_flow;
            stat.explained_flow = decomp_stats.explained_flow;
            stat.weight = decomp_stats.weight;

          
        }
    }
}
fn parse_decomp_file(file_path: &Path) -> std::io::Result<Option<DecompStats>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let mut runtime = 0.0;
    let mut objective_value = 0.0;
    let mut total_flow = 0.0;



    let mut path_weights = Vec::new();
    let mut parsing_paths = false;

for line in reader.lines() {
    let line = line?;

    if line.starts_with("Runtime: ") {
        runtime = line.split_whitespace().nth(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    } else if line.starts_with("Objective Value: ") {
        objective_value = line.split_whitespace().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    } else if line.starts_with("Total Flow: ") {
        total_flow = line.split_whitespace().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    } else if line.starts_with("Paths and Weights:") {
        parsing_paths = true;
    } else if parsing_paths {
        if line.trim().is_empty() {
            parsing_paths = false;
        } else {
            // weight is the first whitespace-separated field
            if let Some(weight_str) = line.split_whitespace().next() {
                if let Ok(weight) = weight_str.parse::<f64>() {
                    path_weights.push(weight);
                }
            }
        }
    }
}

    if runtime > 0.0 || objective_value > 0.0 || total_flow > 0.0 {
        let explained_flow = if total_flow > 0.0 {
            (total_flow - objective_value) / total_flow
        } else {
            0.0
        };
        
        let total_weight: f64 = path_weights.iter().sum();


        Ok(Some(DecompStats {
            runtime,
            objective_value,
            total_flow,
            explained_flow,
            weight: total_weight, 
        }))
    } else {
        Ok(None)
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
    (1, 1)
}

fn parse_alignment_file(
    file_path: &Path, 
    sample_name: String,
    subgraph_name: String,
    part_numbers: (usize, usize)
) -> std::io::Result<AlignmentStats> {
    let file = File::open(file_path)?;
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
        runtime: 0.0,
        objective_value: 0.0,
        nodes: 0,
        edges: 0,
        sources: 0,
        sinks: 0,
        total_flow: 0.0,
        explained_flow: 0.0,    
        weight: 0.0

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

    writer.write_record(&[
        "Sample",
        "Subgraph",
        "Path",
        "Total Paths",
        "Length",
        "Identity %",
        "Identity Count",
        "Gaps %",
        "Gaps Count",
        "Score",
        "Start Position",
        "End Position",
        "Alignment Length",
        "Runtime (s)",
        "Objective Value",
        "Nodes",
        "Edges",
        "Sources (from 0)",
        "Sinks (to 1)",
        "Total Flow",
        "Explained Flow",
        "Weight",
    ])?;

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
            &format!("{:.4}", stats.runtime),
            &format!("{:.6}", stats.objective_value),
            &stats.nodes.to_string(),
            &stats.edges.to_string(),
            &stats.sources.to_string(),
            &stats.sinks.to_string(),
            &format!("{:.6}", stats.total_flow),
            &format!("{:.6}", stats.explained_flow),
            &format!("{:.6}", stats.weight),
        ])?;
    }

    writer.flush()?;
    Ok(())
}
