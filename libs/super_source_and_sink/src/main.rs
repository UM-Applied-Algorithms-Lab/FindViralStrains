use std::collections::HashSet;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

// Helper function to read nodes from a file
fn read_nodes_from_file(filename: &str) -> io::Result<HashSet<String>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let nodes: HashSet<String> = reader
        .lines()
        .filter_map(|line| line.ok())
        .collect();

    Ok(nodes)
}

// Helper function to read edges with weights from a file
fn read_edges_with_weights(filename: &str, output_file: &mut File) -> io::Result<Vec<(String, String, i32)>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut edges = Vec::new();
    let mut line_count = 0;

    
for line in reader.lines() {
    line_count += 1;

    if line_count == 2 {
        // Write the second line (summary line) to the output file
        if let Ok(summary_line) = line { // Ensure line is Ok before parsing
            // Parse the line directly to a u32
            let new_summary_line = summary_line.parse::<u32>().unwrap() + 2; // Add two for the super source and sink nodes
            writeln!(output_file, "# Counts indicate weights on edges")?;
            writeln!(output_file, "{}", new_summary_line)?;
        }
        continue; // Skip further processing for this line
    }
    // Further processing for other lines can go here if needed
        // Process other lines for edges with weights
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 3 {
                let from = parts[0].to_string();
                let to = parts[1].to_string();
                
                // Attempt to parse the weight as an integer
                if let Ok(weight) = parts[2].parse::<i32>() {
                    edges.push((from, to, weight));
                }
            }
        }
    }

    Ok(edges)
}

// Create "super source" and "super sink" nodes and return a list of edges
fn create_super_sources_and_sinks(
    sources_file: &str,
    full_node_list: &str,
    sinks_file: &str,
    edge_file: &str,
    output_file: &mut File,
) -> io::Result<(HashSet<String>, Vec<(String, String, i32)>)> {
    // Read the node lists from the files
    let full_nodes = read_nodes_from_file(full_node_list)?;
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;

    // Read edges with weights and write summary line to output file
    let mut edges = read_edges_with_weights(edge_file, output_file)?;

    // Add edges from the "super source" (node "0") to all source nodes with weight 0
    let super_source = "0".to_string();
    for source in &sources {
        edges.push((super_source.clone(), source.clone(), 0));
    }

    // Add edges from each sink node to the "super sink" (node "1") with weight 0
    let super_sink = "1".to_string();
    for sink in &sinks {
        edges.push((sink.clone(), super_sink.clone(), 0));
    }

    Ok((full_nodes, edges))
}

// Helper function to extract sample name from file path
fn extract_sample_name(file_path: &str) -> String {
    Path::new(file_path)
        .file_stem()
        .and_then(|stem| stem.to_str())
        .unwrap_or("sample")
        .to_string()
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 6 {
        eprintln!("Usage: {} <sources_file> <full_node_list> <sinks_file> <edge_file> <output_directory>", args[0]);
        std::process::exit(1);
    }

    // Command line arguments
    let sources_file = &args[1];
    let full_node_list = &args[2];
    let sinks_file = &args[3];
    let edge_file = &args[4];
    let output_directory = Path::new(&args[5]);

    // Ensure the output directory exists or create it
    if !output_directory.exists() {
        fs::create_dir_all(output_directory)?;
    }

    // Extract sample name from edge file
    let sample_name = extract_sample_name(edge_file);
    let output_file_path = output_directory.join(format!("{}.super.wg", sample_name));

    // Open the output file for writing
    let mut output_file = File::create(&output_file_path)?;

    // Create super sources and sinks, and retrieve the updated nodes and edges list
    let (_full_nodes, edges_with_supers) = create_super_sources_and_sinks(sources_file, full_node_list, sinks_file, edge_file, &mut output_file)?;

    // Append the super nodes and updated edges to the output file
    for (from, to, weight) in &edges_with_supers {
        writeln!(output_file, "{} {} {}", from, to, weight)?;
    }

    println!("New nodes and edges with weights written to: {}", output_file_path.display());
    Ok(())
}
