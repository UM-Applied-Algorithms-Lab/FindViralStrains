use std::collections::HashSet;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

// Helper function to read nodes from a file
fn read_nodes_from_file(filename: &str) -> io::Result<HashSet<String>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let nodes: HashSet<String> = reader.lines().filter_map(|line| line.ok()).collect();

    Ok(nodes)
}

// Helper function to read edges with weights and k-mers from a file
fn read_edges_with_weights(
    filename: &str,
    output_file: &mut File,
) -> io::Result<(HashSet<String>, Vec<(String, String, i32, String)>)> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut edges = Vec::new();
    let mut nodes = HashSet::new();
    let mut line_count = 0;

    for line in reader.lines() {
        line_count += 1;

        if line_count == 2 {
            // Write the second line (summary line) to the output file
            if let Ok(summary_line) = line {
                // Ensure line is Ok before parsing
                // Parse the line directly to a u32
                let new_summary_line = summary_line.parse::<u32>().unwrap() + 2; // Add two for the super source and sink nodes
                writeln!(output_file, "# Counts indicate weights on edges")?;
                writeln!(output_file, "{}", new_summary_line)?;
            }
            continue;
        }
        // Further processing for other lines can go here if needed
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 4 {
                let from = parts[0].to_string();
                let to = parts[1].to_string();

                // Add nodes to the full node list
                nodes.insert(from.clone());
                nodes.insert(to.clone());

                // Attempt to parse the weight as an integer
                if let Ok(weight) = parts[2].parse::<i32>() {
                    let kmer = parts[3].to_string();
                    edges.push((from, to, weight, kmer));
                }
            }
        }
    }

    Ok((nodes, edges))
}

// Create "super source" and "super sink" nodes and return a list of edges
fn create_super_sources_and_sinks(
    sources_file: &str,
    sinks_file: &str,
    edge_file: &str,
    output_file: &mut File,
) -> io::Result<(HashSet<String>, Vec<(String, String, i32, String)>)> {
    // Read the node lists from the files
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;

    // Read edges and derive full nodes from the edge file
    let (full_nodes, mut edges) = read_edges_with_weights(edge_file, output_file)?;

    // Add edges from the "super source" (node "0") to all source nodes with weight 0
    let super_source = "0".to_string();
    for source in &sources {
        edges.push((super_source.clone(), source.clone(), 0, "".to_string()));
    }

    // Add edges from each sink node to the "super sink" (node "1") with weight 0
    let super_sink = "1".to_string();
    for sink in &sinks {
        edges.push((sink.clone(), super_sink.clone(), 0, "".to_string()));
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

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 5 {
        eprintln!(
            "Usage: {} <sources_file> <sinks_file> <edge_file> <output_directory>",
            args[0]
        );
        std::process::exit(1);
    }

    // Command line arguments
    let sources_file = &args[1];
    let sinks_file = &args[2];
    let edge_file = &args[3];
    let output_directory = Path::new(&args[4]);

    // Ensure the output directory exists or create it
    if !output_directory.exists() {
        fs::create_dir_all(output_directory).expect("unable to generate out directories");
    }

    // Extract sample name from edge file
    let sample_name = extract_sample_name(edge_file);
    let output_file_path = output_directory.join(format!("{}.super.wg", sample_name));

    // Open the output file for writing
    let mut output_file = File::create(&output_file_path).expect(&format!(
        "unable to create output file at location {}",
        output_file_path.display()
    ));

    // Create super sources and sinks, and retrieve the updated nodes and edges list
    let (_full_nodes, edges_with_supers) = create_super_sources_and_sinks(
        sources_file,
        sinks_file,
        edge_file,
        &mut output_file,
    )
    .expect("unable to create super sources and sinks");

    // Append the super nodes and updated edges to the output file
    for (from, to, weight, kmer) in &edges_with_supers {
        writeln!(output_file, "{} {} {} {}", from, to, weight, kmer).expect(&format!(
            "failure in writing edge data to output file at location {}",
            output_file_path.display()
        ));
    }

    println!(
        "New nodes and edges with weights written to: {}",
        output_file_path.display()
    );
}