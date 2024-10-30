use std::collections::HashSet;
use std::env;
use std::io::{self, BufRead, BufReader, Write};
use std::fs::{File};

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
) -> io::Result<(HashSet<String>, Vec<(String, String)>)> {
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

    Ok((full_nodes, edges))
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

    // Create super sources and sinks and retrieve the updated nodes and edges list
    let (full_nodes, edges) = create_super_sources_and_sinks(sources_file, full_node_list, sinks_file)?;

    // Create a new file name based on the full_node_list file name with ".super.mg" suffix
    let mut base_name = full_node_list.to_string();
    if base_name.ends_with(".mg") {
        // Remove the ".mg" suffix if it exists
        base_name = base_name.trim_end_matches(".mg").to_string();
    }

    let output_file_name = format!("{}.super.mg", base_name);
    let mut output_file = File::create(&output_file_name)?;

    // Print all nodes including full nodes and super nodes
    writeln!(output_file, "# Full node list:")?;
    for node in &full_nodes {
        writeln!(output_file, "{}", node)?;
    }
    writeln!(output_file, "0")?; // Super Source Node
    writeln!(output_file, "1")?; // Super Sink Node
    writeln!(output_file)?;

    // Print and append the edges to the new file
    writeln!(output_file, "# New edges created:")?;
    for (from, to) in &edges {
        writeln!(output_file, "{} {}", from, to)?;
    }

    println!("New nodes and edges written to: {}", output_file_name);
    Ok(())
}
