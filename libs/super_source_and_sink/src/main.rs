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
) -> io::Result<(HashSet<String>, Vec<(String, String, i32, String)>)> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut edges = Vec::new();
    let mut nodes = HashSet::new();
    let mut line_count = 0;

    for line in reader.lines() {
        line_count += 1;

        // Skip the first line (header) and process the rest
        if line_count == 1 {
            continue;
        }

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
) -> io::Result<()> {
    // Read the node lists from the files
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;

    // Read edges and derive full nodes from the edge file
    let (mut full_nodes, edges) = read_edges_with_weights(edge_file)?;

    // Add super source ("0") and super sink ("1") to the full node list
    full_nodes.insert("0".to_string());
    full_nodes.insert("1".to_string());

    // Write all original edges to the output file
    for (from, to, weight, kmer) in &edges {
        writeln!(output_file, "{} {} {} {}", from, to, weight, kmer)?;
    }

    // Add edges from the "super source" (node "0") to all source nodes with weight 0
    let super_source = "0".to_string();
    for source in &sources {
        writeln!(output_file, "{} {} 0", super_source, source)?;
    }

    // Add edges from each sink node to the "super sink" (node "1") with weight 0
    let super_sink = "1".to_string();
    for sink in &sinks {
        writeln!(output_file, "{} {} 0", sink, super_sink)?;
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 5 {
        eprintln!(
            "Usage: {} <sources_file> <sinks_file> <edge_file> <output_file>",
            args[0]
        );
        std::process::exit(1);
    }

    // Command line arguments
    let sources_file = &args[1];
    let sinks_file = &args[2];
    let edge_file = &args[3];
    let output_file_path = Path::new(&args[4]);

    // Ensure the output directory exists or create it
    if let Some(parent) = output_file_path.parent() {
        if !parent.exists() {
            fs::create_dir_all(parent).expect("unable to generate output directories");
        }
    }

    // Open the output file for writing
    let mut output_file = File::create(&output_file_path).expect(&format!(
        "unable to create output file at location {}",
        output_file_path.display()
    ));

    // Copy the entire content of the edge file to the output file
    let edge_file_content = fs::read_to_string(edge_file).expect("unable to read edge file");
    write!(output_file, "{}", edge_file_content).expect("unable to write edge file content");

    // Create super sources and sinks, and write all edges to the output file
    create_super_sources_and_sinks(
        sources_file,
        sinks_file,
        edge_file,
        &mut output_file,
    )
    .expect("unable to create super sources and sinks");

    println!(
        "New nodes and edges with weights written to: {}",
        output_file_path.display()
    );
}