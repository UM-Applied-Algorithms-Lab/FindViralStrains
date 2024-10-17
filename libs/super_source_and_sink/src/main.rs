use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

fn read_nodes_from_file(filename: &str) -> io::Result<HashSet<String>> {
    let path = Path::new(filename);
    let file = File::open(&path)?;
    let reader = io::BufReader::new(file);
    
    let nodes: HashSet<String> = reader
        .lines()
        .filter_map(Result::ok)
        .collect();
        
    Ok(nodes)
}

fn create_super_sources_and_sinks(
    sources_file: &str,
    full_node_list: &str,
    sinks_file: &str,
) -> io::Result<(HashSet<String>, HashSet<String>)> {
    // Read the node lists from the files
    let full_nodes = read_nodes_from_file(full_node_list)?;
    let sinks = read_nodes_from_file(sinks_file)?;
    let sources = read_nodes_from_file(sources_file)?;
    
    // Create super source and super sink sets
    let mut super_source = HashSet::new();
    let mut super_sink = HashSet::new();
    
    // Super source: nodes that are sources but not sinks
    for node in sources.difference(&sinks) {
        if full_nodes.contains(node) {
            super_source.insert(node.clone());
        }
    }
    
    // Super sink: nodes that are sinks but not sources
    for node in sinks.difference(&sources) {
        if full_nodes.contains(node) {
            super_sink.insert(node.clone());
        }
    }

    Ok((super_source, super_sink))
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        eprintln!("Usage: {} <sources_file> <full_node_list> <sinks_file>", args[0]);
        std::process::exit(1);
    }

    let sources_file = &args[1];
    let full_node_list = &args[2];
    let sinks_file = &args[3];

    let (super_source, super_sink) = create_super_sources_and_sinks(sources_file, full_node_list, sinks_file)?;
    
    println!("Super Sources: {:?}", super_source);
    println!("Super Sinks: {:?}", super_sink);
    
    Ok(())
}
