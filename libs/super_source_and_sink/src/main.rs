use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::hash::{BuildHasher, Hasher};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

// Custom hasher with fixed seed
#[derive(Default)]
struct FixedHasher(u64);

impl Hasher for FixedHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        // Simple deterministic hash function - DJB2 algorithm
        let mut hash: u64 = 5381;
        for &byte in bytes {
            hash = ((hash << 5).wrapping_add(hash)).wrapping_add(byte as u64);
        }
        self.0 = hash;
    }
}

#[derive(Default)]
struct FixedBuildHasher;

impl BuildHasher for FixedBuildHasher {
    type Hasher = FixedHasher;

    fn build_hasher(&self) -> FixedHasher {
        FixedHasher(123) // Fixed seed
    }
}

// Helper function to read nodes from a file
fn read_nodes_from_file(filename: &str) -> io::Result<HashMap<String, (), FixedBuildHasher>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut nodes = HashMap::with_hasher(FixedBuildHasher);
    for line in reader.lines() {
        if let Ok(node) = line {
            nodes.insert(node, ());
        }
    }

    Ok(nodes)
}

// Helper function to read edges with weights and k-mers from a file
fn read_edges_with_weights(
    filename: &str,
) -> io::Result<(
    HashMap<String, (), FixedBuildHasher>,
    Vec<(String, String, i32)>,
)> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut edges = Vec::new();
    let mut nodes = HashMap::with_hasher(FixedBuildHasher);
    let mut line_count = 0;

    for line in reader.lines() {
        line_count += 1;

        // Skip the first line (header) and process the rest
        if line_count == 1 {
            continue;
        }

        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 3 {
                let from = parts[0].to_string();
                let to = parts[1].to_string();

                // Add nodes to the full node list
                nodes.insert(from.clone(), ());
                nodes.insert(to.clone(), ());

                // Attempt to parse the weight as an integer
                if let Ok(weight) = parts[2].parse::<i32>() {
                    edges.push((from, to, weight));
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
    graph_name: &str,
) -> io::Result<()> {
    // Write the graph name as a comment at the top of the output file
    writeln!(output_file, "# {}", graph_name)?;

    // Read the node lists from the files
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;

    // Read edges and derive full nodes from the edge file
    let (mut full_nodes, edges) = read_edges_with_weights(edge_file)?;

    // Add super source ("0") and super sink ("1") to the full node list
    full_nodes.insert("0".to_string(), ());
    full_nodes.insert("1".to_string(), ());

    // Read and sort all original edges for deterministic output
    let edge_file_content = fs::read_to_string(edge_file).expect("unable to read edge file");
    let mut edge_lines: Vec<&str> = edge_file_content.lines().collect();

    // Skip the first line (header) and sort the remaining edges
    if !edge_lines.is_empty() {
        let header = edge_lines.remove(0);
        edge_lines.sort(); // Sort edges for deterministic output
        edge_lines.insert(0, header); // Put header back at the beginning
    }

    // Write all sorted original edges to the output file
    for line in edge_lines {
        writeln!(output_file, "{}", line)?;
    }

    // Convert HashMap keys to sorted vectors for deterministic output
    let mut sorted_sources: Vec<String> = sources.keys().cloned().collect();
    sorted_sources.sort();

    let mut sorted_sinks: Vec<String> = sinks.keys().cloned().collect();
    sorted_sinks.sort();

    // Add edges from the "super source" (node "0") to all source nodes with weight 0
    let super_source = "0".to_string();
    for source in sorted_sources {
        writeln!(output_file, "{} {} 0", super_source, source)?;
    }

    // Add edges from each sink node to the "super sink" (node "1") with weight 0
    let super_sink = "1".to_string();
    for sink in sorted_sinks {
        writeln!(output_file, "{} {} 0", sink, super_sink)?;
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 6 {
        eprintln!(
            "Usage: {} <sources_file> <sinks_file> <edge_file> <output_file> <graph_name>",
            args[0]
        );
        std::process::exit(1);
    }

    // Command line arguments
    let sources_file = &args[1];
    let sinks_file = &args[2];
    let edge_file = &args[3];
    let output_file_path = Path::new(&args[4]);
    let graph_name = &args[5];

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

    // Create super sources and sinks, and write all edges to the output file
    create_super_sources_and_sinks(
        sources_file,
        sinks_file,
        edge_file,
        &mut output_file,
        graph_name,
    )
    .expect("unable to create super sources and sinks");
}
