use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Result};
use std::path::Path;

// Structure to hold graph analysis data //
struct GraphAnalysisData {
    num_edges: usize,                 // Total number of edges in the graph //
    num_disconnected_subgraphs: usize, // Number of disconnected subgraphs //
    sources: Vec<String>,             // List of source nodes (nodes with no incoming edges) //
    sinks: Vec<String>,               // List of sink nodes (nodes with no outgoing edges) //
}

fn main() {
    // Collects command-line arguments into a vector //
    let program_args: Vec<String> = env::args().collect();
    
    // Checks if there are two arguments (program name and file path) //
    if env::args().len() < 2 {
        println!("Error: requires at least 1 arg: file src for .mg file to analyze.");
    } else {
        let mg_file_src: &str = &program_args[1];

        // Analyzes the graph and prints out the results if successful //
        match make_graph_analysis_data(Path::new(mg_file_src)) {
            Ok(stats) => println!(
                "   num edges,\t  num source nodes,\t num sink nodes,\tnum disconnected graphs\n {}\t,{}\t,{},\t{}",
                stats.num_edges, stats.sources.len(), stats.sinks.len(), stats.num_disconnected_subgraphs),
            Err(e) => println!("could not process input file: {}", e.to_string()),
        };
    }
}

// Processes the file and returns graph analysis data //
fn make_graph_analysis_data(file_path: &Path) -> Result<GraphAnalysisData> {
    let edge_map = make_edge_map(file_path)?;
    let (sources, sinks) = make_source_sink_lists(&edge_map);
    let num_disconnected_subgraphs = get_num_disconnected_subgraphs(&edge_map, &sources)?;

    // Calculates the total number of edges //
    let num_edges: usize = edge_map
        .iter()
        .map(|element| element.1.len())
        .sum::<usize>();

    Ok(GraphAnalysisData {
        num_edges,
        num_disconnected_subgraphs,
        sources,
        sinks,
    })
}

// Reads file and creates a map of edges between nodes //
fn make_edge_map(file_path: &Path) -> Result<HashMap<String, Vec<String>>> {
    if !file_path.exists() {
        println!("file does not exist!");
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "File does not exist",
        ));
    }
    let file = File::open(file_path)?;
    let file_reader = BufReader::new(file);

    let mut edge_map: HashMap<String, Vec<String>> = HashMap::new();
    let lines = file_reader.lines().skip(2); // Skips the first two lines, assuming header

    for line in lines {
        match line {
            Ok(line) => {
                let mut split_line = line.split_whitespace();
                // Extracts the "from" node and "to" node from each line
                let from_node = match split_line.next() {
                    Some(node) => node.to_string(),
                    None => {
                        println!("couldn't parse line {}", line);
                        continue;
                    }
                };
                let to_node = match split_line.next() {
                    Some(node) => node.to_string(),
                    None => {
                        println!("couldn't parse line {}", line);
                        continue;
                    }
                };
                // Adds the edge to the map
                edge_map
                    .entry(from_node)
                    .or_insert_with(Vec::new)
                    .push(to_node);
            }
            Err(_) => println!("unable to get file line."),
        }
    }

    Ok(edge_map)
}

// Finds the number of disconnected subgraphs in the graph
fn get_num_disconnected_subgraphs(
    edge_map: &HashMap<String, Vec<String>>,
    sources: &Vec<String>,
) -> Result<usize> {
    let mut encountered_sink_set: HashSet<String> = HashSet::new();
    let mut encountered_nodes: HashSet<String> = HashSet::new();

    let mut num_disconnected_subgraphs = 0;

    for source in sources {
        let mut this_sources_sinks: Vec<String> = Vec::new();
        let mut node_stack: Vec<String> = Vec::new();

        node_stack.push(source.clone());

        while !node_stack.is_empty() {
            let current_node = node_stack.pop().unwrap();

            if !encountered_nodes.contains(&current_node) {
                match edge_map.get(&current_node) {
                    Some(child_nodes) => {
                        for child_node in child_nodes {
                            node_stack.push(child_node.clone());
                        }
                        encountered_nodes.insert(current_node);
                    }
                    None => {
                        // If the node has no outgoing edges, it's a sink
                        this_sources_sinks.push(current_node);
                    }
                }
            }
        }

        // Checks if the sinks from this source are new, indicating a new disconnected subgraph
        let new_sinks: Vec<String> = this_sources_sinks
            .iter()
            .filter(|node| !encountered_sink_set.contains(*node))
            .cloned()
            .collect();

        if new_sinks.len() != 0 {
            num_disconnected_subgraphs += 1;
            for sink in new_sinks {
                encountered_sink_set.insert(sink);
            }
        }
    }

    Ok(num_disconnected_subgraphs)
}

// Identifies source nodes (no incoming edges) and sink nodes (no outgoing edges)
fn make_source_sink_lists(edge_map: &HashMap<String, Vec<String>>) -> (Vec<String>, Vec<String>) {
    let end_nodes: HashSet<&String> = edge_map.values().flatten().collect();

    let sources: Vec<String> = edge_map
        .iter()
        .filter(|node| !end_nodes.contains(node.0))
        .map(|node| node.0)
        .cloned()
        .collect();

    let sinks = end_nodes
        .iter()
        .filter(|node| !edge_map.contains_key(**node))
        .map(|node| *node)
        .cloned()
        .collect();

    (sources, sinks)
}
