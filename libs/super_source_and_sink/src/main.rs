use std::collections::HashSet;
use std::env;
use std::io::{self, BufRead, BufReader, Write};
use std::fs::{File, OpenOptions};

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

// Create "super source" and "super sink" nodes and return a list of edges
fn create_super_sources_and_sinks(
    sources_file: &str,
    full_node_list: &str,
    sinks_file: &str,
) -> io::Result<Vec<(String, String)>> {
    // Read the node lists from the files
    let full_nodes = read_nodes_from_file(full_node_list)?;
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;

    let mut edges = Vec::new();

    // Add edges from the "super source" (node "0") to all source nodes
    let super_source = "0".to_string();
    for source in &sources {
        edges.push((super_source.clone(), source.clone()));
    }

    // Add edges from each sink node to the "super sink" (node "1")
    let super_sink = "1".to_string();
    for sink in &sinks {
        edges.push((sink.clone(), super_sink.clone()));
    }

    Ok(edges)
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        eprintln!("Usage: {} <sources_file> <full_node_list> <sinks_file>", args[0]);
        std::process::exit(1);
    }

    // Command line arguments
    let sources_file = &args[1];
    let full_node_list = &args[2];
    let sinks_file = &args[3];

    // Create super sources and sinks and retrieve the updated edges list
    let edges = create_super_sources_and_sinks(sources_file, full_node_list, sinks_file)?;

    // Print information about the new nodes (super source and super sink)
    println!("New nodes added:");
    println!("  Super Source Node: 0");
    println!("  Super Sink Node: 1");
    println!(); // Line break for better readability

    // Open the full_node_list file in append mode
    let mut file = OpenOptions::new()
        .append(true)
        .open(full_node_list)?;

    // Print the edges in the desired format and append them to the file
    println!("New edges created:");
    for (index, (from, to)) in edges.iter().enumerate() {
        // Example format for printing:
        // index + 1, from_node, to_node, nucleotide_sequence (empty in this case)
        println!("{:>3} {}", from, to);

        // Append each 'from' and 'to' section to the full_node_list file
        writeln!(file, "{} {}", from, to)?;
    }

    Ok(())
}