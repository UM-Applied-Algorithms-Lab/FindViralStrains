use clap::Parser;
use colored::Colorize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::{BuildHasher, Hasher};
use std::io::{BufRead, BufReader, Result, Write};
use std::path::Path;
use std::rc::Rc;

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

//Struct used to handle input args, with the Clap rust crate.
//This pulls help text from the comments, and compiler flags/options from the variable names
#[derive(Parser, Debug)]
struct InputArgs {
    ///file src for the dbg graph file to parse
    #[arg(short, long)]
    dbg_file_name: String,

    ///percentage of nodes required in a subgraph to be reported, e.g., 0.05 requires 5% of total nodes in a subgraph to be reported.
    #[arg(short, long, default_value_t = 0.01)]
    node_percent_cutoff: f32,

    ///displays all subgraphs, not just those that meet the cutoff node percentage
    #[arg(short, long, default_value_t = false)]
    all_subgraphs_displayed: bool,

    /// the directory that the output will be written to. By default, outputs will be put in the calling_directory/input_name/
    #[arg(short, long)]
    output_directory: Option<String>,

    /// file to save statistics output (in addition to terminal display)
    #[arg(short, long)]
    stats_output_file: Option<String>,

    /// exclude graphs with cycles from output
    #[arg(short = 'x', long, default_value_t = false)]
    exclude_cyclic_graphs: bool,
}

//struct to hold the analysis data for a graph or subgraph
struct GraphAnalysisData {
    num_nodes: usize,
    num_edges: usize,
    total_edge_weight: usize,
    num_disconnected_subgraphs: usize,
    sources: Vec<Rc<str>>,
    sinks: Vec<Rc<str>>,
    is_acyclic: bool,
}

enum SubgraphDisplayType {
    All,
    OnlySignificant(f32),
}

impl std::fmt::Display for GraphAnalysisData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:<12}{:<12}{:<18}{:<14}{:<12}{:<12}{:<10}\n{:<12}{:<12}{:<18}{:<14}{:<12}{:<12}{:<10}",
            "nodes".blue(),
            "edges".blue(),
            "total_weight".blue(),
            "subgraphs".blue(),
            "sources".blue(),
            "sinks".blue(),
            "acyclic?".blue(),
            self.num_nodes,
            self.num_edges,
            self.total_edge_weight,
            self.num_disconnected_subgraphs,
            self.sources.len(),
            self.sinks.len(),
            match self.is_acyclic {
                true => "yes".green(),
                false => "no".red(),
            }
        )
    }
}

//for any node in a graph, store the in and out edges as vectors
struct NodeEdges {
    in_edges: Vec<Rc<str>>,
    out_edges: Vec<Rc<str>>,
}

impl Clone for NodeEdges {
    fn clone(&self) -> Self {
        Self {
            in_edges: self.in_edges.clone(),
            out_edges: self.out_edges.clone(),
        }
    }
}

impl NodeEdges {
    pub fn new() -> Self {
        Self {
            in_edges: Vec::new(),
            out_edges: Vec::new(),
        }
    }
}

fn main() {
    let args = InputArgs::parse();

    //parse the main de bruijn graph from the input mer-graph file
    let (main_graph, edge_info, graph_label) = match make_main_graph(Path::new(&args.dbg_file_name))
    {
        Ok(graph) => graph,
        Err(err) => panic!("Unable to generate graph, check file: {}", err),
    };

    //generate the list of subgraphs
    let subgraph_list = make_subgraph_list(&main_graph);
    let significant_subgraph_list =
        make_significant_subgraph_list(&subgraph_list, args.node_percent_cutoff);

    //displays the stats from the main graph and any subgraphs that match the display criteria
    let display_type = match args.all_subgraphs_displayed {
        true => SubgraphDisplayType::All,
        false => SubgraphDisplayType::OnlySignificant(args.node_percent_cutoff),
    };
    if let Err(e) = display_graph_stats(
        &main_graph,
        &subgraph_list,
        &edge_info,
        display_type,
        &args.stats_output_file,
    ) {
        eprintln!("Warning: Failed to write statistics to file: {}", e);
    }

    //writes the subgraphs over the node % cutoff to individual files for further processing
    write_subgraph_files(
        significant_subgraph_list,
        &args.dbg_file_name,
        &graph_label,
        &edge_info,
        args.output_directory,
        args.exclude_cyclic_graphs,
    );
}

/// given a list of subgraphs, writes them to files for further use in the pipeline.
/// writes subgraph source, sink, and .mg files to a subdirectory named after the input .mg file src
///
/// significant_subgraph_list should be the list of all subgraphs to print
fn write_subgraph_files(
    significant_subgraph_list: Vec<&HashMap<Rc<str>, NodeEdges, FixedBuildHasher>>,
    base_file_name: &String,
    main_graph_label: &String,
    edge_info: &HashMap<(Rc<str>, Rc<str>), (usize, Rc<str>), FixedBuildHasher>,
    output_dir: Option<String>,
    exclude_cyclic_graphs: bool,
) {
    for (subgraph_idx, subgraph) in significant_subgraph_list.iter().enumerate() {
        // Skip cyclic graphs if the flag is set
        if exclude_cyclic_graphs && !graph_is_acyclic(subgraph) {
            continue;
        }

        let subgraph_sub_dir = base_file_name.to_string() + "_subgraphs";
        let subgraph_directory_name = match output_dir {
            Some(ref dir) => Path::new(dir).join(Path::new(&subgraph_sub_dir)),
            None => Path::new(&subgraph_sub_dir).to_path_buf(),
        };

        std::fs::create_dir_all(&subgraph_directory_name)
            .expect("could not create subgraph directory");
        let subgraph_idx_string = subgraph_idx.to_string();

        let (mut sources, mut sinks) = make_source_sink_lists(&subgraph);
        // Sort sources and sinks for deterministic output
        sources.sort();
        sinks.sort();

        let _ = std::fs::write(
            Path::new(&subgraph_directory_name)
                .join(format!("graph_{}.sources", &subgraph_idx_string)),
            sources.join("\n"),
        );
        let _ = std::fs::write(
            Path::new(&subgraph_directory_name)
                .join(format!("graph_{}.sinks", &subgraph_idx_string)),
            sinks.join("\n"),
        );

        let subgraph_file_path =
            Path::new(&subgraph_directory_name).join(format!("graph_{}.dbg", &subgraph_idx_string));

        //configure the file output to initially overwrite the file, and append all nodes of the subgraph
        let mut subgraph_mg_file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true) // Overwrite if exists
            .open(subgraph_file_path)
            .expect("unable to open subgraph mg file for writing");

        //write the subgraph label
        subgraph_mg_file
            .write_fmt(format_args!(
                "{} subgraph_{}\n",
                main_graph_label, subgraph_idx
            ))
            .expect("unable to write subgraph label to file");

        //get the number of nodes for the subgraph
        let subgraph_num_nodes = subgraph.len();
        subgraph_mg_file
            .write_fmt(format_args!("{}\n", subgraph_num_nodes))
            .expect("unable to write subgraph node count to file");

        // Sort nodes for deterministic output
        let mut sorted_nodes: Vec<&Rc<str>> = subgraph.keys().collect();
        sorted_nodes.sort();

        for from_node in sorted_nodes {
            let edges = subgraph.get(from_node).unwrap();
            // Sort out_edges for deterministic output
            let mut sorted_out_edges: Vec<&Rc<str>> = edges.out_edges.iter().collect();
            sorted_out_edges.sort();

            for to_node in sorted_out_edges {
                match edge_info.get(&(from_node.clone(), to_node.clone())) {
                    Some(result) => {
                        let (count, kmer) = result;
                        subgraph_mg_file
                            .write_fmt(format_args!(
                                "{}\t{}\t{}\t{}\n",
                                from_node, to_node, count, kmer
                            ))
                            .expect("unable to write graph line to subgraph file");
                    }
                    None => {
                        println!("No count or kmer stored ")
                    }
                }
            }
        }
    }
}

/// displays the statistics of the main graph and all subgraphs (or only all significant subgraphs)
fn display_graph_stats(
    main_graph: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
    subgraph_list: &Vec<HashMap<Rc<str>, NodeEdges, FixedBuildHasher>>,
    edge_info: &HashMap<(Rc<str>, Rc<str>), (usize, Rc<str>), FixedBuildHasher>,
    subgraph_display_type: SubgraphDisplayType,
    stats_output_file: &Option<String>,
) -> std::io::Result<()> {
    let mut output = String::new();

    output.push_str(&format!(
        "Main Graph Stats:\n{}\n",
        make_graph_stats(main_graph, &edge_info, subgraph_list.len())
    ));

    match subgraph_display_type {
        SubgraphDisplayType::All => {
            // Create indexed list and sort by size then original index for deterministic order
            let mut indexed_subgraphs: Vec<(
                usize,
                &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
            )> = subgraph_list.iter().enumerate().collect();

            // Sort by size (descending) then by original index (ascending) for deterministic order
            indexed_subgraphs.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then(a.0.cmp(&b.0)));

            for (original_idx, subgraph) in indexed_subgraphs {
                output.push_str(&format!(
                    "Subgraph {}:\n{}\n",
                    original_idx,
                    make_graph_stats(subgraph, &edge_info, 0)
                ));
            }
        }

        SubgraphDisplayType::OnlySignificant(cutoff) => {
            let significant_subgraph_list = make_significant_subgraph_list(&subgraph_list, cutoff);
            for (subgraph_idx, subgraph) in significant_subgraph_list.iter().enumerate() {
                output.push_str(&format!(
                    "Subgraph {}:\n{}\n",
                    subgraph_idx,
                    make_graph_stats(subgraph, &edge_info, 0)
                ));
            }
        }
    }

    // Print to terminal (original functionality)
    print!("{}", output);

    // Write to file if specified
    if let Some(file_path) = stats_output_file {
        std::fs::write(file_path, output)?;
    }

    Ok(())
}

/// generates the stats for a graph (or subgraph)
/// If using for the full graph, give the num_subgraphs, for a subgraph just give zero, I guess...
fn make_graph_stats(
    graph: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
    edge_info: &HashMap<(Rc<str>, Rc<str>), (usize, Rc<str>), FixedBuildHasher>,
    num_subgraphs: usize,
) -> GraphAnalysisData {
    let (mut sources, mut sinks) = make_source_sink_lists(&graph);
    // Sort sources and sinks for deterministic output
    sources.sort();
    sinks.sort();

    // Calculate total edge weight with deterministic iteration order
    let mut sorted_nodes: Vec<&Rc<str>> = graph.keys().collect();
    sorted_nodes.sort();

    let mut total_edge_weight = 0;
    for from_node in &sorted_nodes {
        let edges = graph.get(*from_node).unwrap();
        let mut sorted_out_edges: Vec<&Rc<str>> = edges.out_edges.iter().collect();
        sorted_out_edges.sort();

        for to_node in &sorted_out_edges {
            if let Some((count, _)) = edge_info.get(&((*from_node).clone(), (*to_node).clone())) {
                total_edge_weight += count;
            }
        }
    }

    // Calculate number of edges with deterministic iteration
    let num_edges: usize = sorted_nodes
        .iter()
        .map(|node| graph.get(*node).unwrap().out_edges.len())
        .sum();

    GraphAnalysisData {
        num_nodes: graph.len(),
        num_edges,
        total_edge_weight,
        num_disconnected_subgraphs: num_subgraphs,
        sources,
        sinks,
        is_acyclic: graph_is_acyclic(&graph),
    }
}

/// generates the list of subgraphs. This is a new graph list, so it increaes memory, but avoids memory issues and the
/// borrow checker (and lifetimes)
fn make_subgraph_list(
    main_graph: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
) -> Vec<HashMap<Rc<str>, NodeEdges, FixedBuildHasher>> {
    let mut node_colors: HashMap<Rc<str>, usize, FixedBuildHasher> =
        HashMap::with_hasher(FixedBuildHasher);
    let mut subgraph_list: Vec<HashMap<Rc<str>, NodeEdges, FixedBuildHasher>> = Vec::new();

    // Sort nodes for deterministic flood-fill order
    let mut sorted_nodes: Vec<&Rc<str>> = main_graph.keys().collect();
    sorted_nodes.sort();

    // flood fills the graph to find all connected nodes
    for (node_idx, node_name) in sorted_nodes.iter().enumerate() {
        graph_color_flood_fill(main_graph, &mut node_colors, node_name, node_idx);
    }

    // finds the number of colors that were used to color all nodes of the graph, and therefore,
    // finds the set (and count) of subgraphs
    let mut node_color_set: Vec<usize> = node_colors
        .iter()
        .map(|(_, node_color)| *node_color)
        .collect();
    node_color_set.sort(); // Sort colors for deterministic subgraph order
    node_color_set.dedup();

    // generates separate hashmaps for each subgraph
    for node_color in node_color_set {
        let mut subgraph_entries: Vec<(Rc<str>, NodeEdges)> = node_colors
            .iter()
            .filter(|(_, color)| node_color == **color)
            .map(|(node_name, _)| {
                (
                    node_name.clone(),
                    main_graph.get(node_name).unwrap().clone(),
                )
            })
            .collect();

        // Sort subgraph entries for deterministic HashMap creation
        subgraph_entries.sort_by(|a, b| a.0.cmp(&b.0));

        let subgraph: HashMap<Rc<str>, NodeEdges, FixedBuildHasher> =
            subgraph_entries.into_iter().collect();
        subgraph_list.push(subgraph);
    }

    // Sort subgraph list by size (largest first) then by minimum node name for consistent ordering
    subgraph_list.sort_by(|a, b| {
        b.len().cmp(&a.len()).then_with(|| {
            let a_min = a.keys().min().unwrap();
            let b_min = b.keys().min().unwrap();
            a_min.cmp(b_min)
        })
    });

    return subgraph_list;
}

/// using a percent cutoff, filters the list of subgraphs to only those who have a large enough share of nodes.
/// "large enough share" is defined as at least 'node_percent_cutoff' percent of the total nodes of the full graph.
fn make_significant_subgraph_list(
    subgraph_list: &Vec<HashMap<Rc<str>, NodeEdges, FixedBuildHasher>>,
    node_percent_cutoff: f32,
) -> Vec<&HashMap<Rc<str>, NodeEdges, FixedBuildHasher>> {
    let num_total_nodes: usize = subgraph_list.iter().map(|subgraph| subgraph.len()).sum();
    let min_nodes_for_significant_subgraph: usize =
        ((num_total_nodes as f32) * node_percent_cutoff) as usize;

    let mut significant: Vec<&HashMap<Rc<str>, NodeEdges, FixedBuildHasher>> = subgraph_list
        .iter()
        .filter(|subgraph| subgraph.len() >= min_nodes_for_significant_subgraph)
        .collect();

    // Sort significant subgraphs by size (descending) then by minimum node name for deterministic order
    significant.sort_by(|a, b| {
        b.len().cmp(&a.len()).then_with(|| {
            let a_min = a.keys().min().unwrap();
            let b_min = b.keys().min().unwrap();
            a_min.cmp(b_min)
        })
    });

    return significant;
}

/// reads in the mer-graph from the file_path, and generates the full de bruijn graph as a hashmap
/// of outgoing edges and a hashmap of the kmers labeling the edges
fn make_main_graph(
    file_path: &Path,
) -> Result<(
    HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
    HashMap<(Rc<str>, Rc<str>), (usize, Rc<str>), FixedBuildHasher>,
    String,
)> {
    let file = File::open(file_path)?;
    let file_reader = BufReader::new(file);

    let mut main_graph: HashMap<Rc<str>, NodeEdges, FixedBuildHasher> =
        HashMap::with_hasher(FixedBuildHasher);
    let mut edge_info: HashMap<(Rc<str>, Rc<str>), (usize, Rc<str>), FixedBuildHasher> =
        HashMap::with_hasher(FixedBuildHasher);
    let lines = file_reader.lines();
    // decompose needs a graph label, but assembly graph generator does not make one //
    let graph_label = "# fake label".to_string();

    for line in lines {
        let line = line.expect("unable to read line in input file");
        let mut split_line = line.split_whitespace();

        let from_node: Rc<str> = Rc::from(split_line.next().expect(
            "encountered incorrectly formatted line in input file: \
            could not parse first node name",
        ));

        let to_node: Rc<str> = Rc::from(split_line.next().expect(
            "encountered incorrectly formatted line in input file: \
            could not parse second node name",
        ));

        let count: usize = split_line.next().unwrap().parse().expect(
            // TODO change to int
            "encountered incorrectly formatted line in input file: \
            could not parse second node name",
        );

        let edge_kmer: Rc<str> = Rc::from(split_line.next().expect(
            "encountered incorrectly formatted line in input file: \
            could not parse second node name",
        ));

        main_graph
            .entry(from_node.clone())
            .or_insert_with(NodeEdges::new)
            .out_edges
            .push(to_node.clone());
        main_graph
            .entry(to_node.clone())
            .or_insert_with(NodeEdges::new)
            .in_edges
            .push(from_node.clone());

        edge_info.insert((from_node, to_node), (count, edge_kmer));
    }

    return Ok((main_graph, edge_info, graph_label));
}

/// generates lists of sources and sinks for the given graph or subgraph
fn make_source_sink_lists(
    edge_map: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
) -> (Vec<Rc<str>>, Vec<Rc<str>>) {
    let sources: Vec<Rc<str>> = edge_map
        .iter()
        .filter(|(_, edges)| edges.in_edges.is_empty())
        .map(|(node_name, _)| node_name.clone())
        .collect();

    let sinks: Vec<Rc<str>> = edge_map
        .iter()
        .filter(|(_, edges)| edges.out_edges.is_empty())
        .map(|(node_name, _)| node_name.clone())
        .collect();

    (sources, sinks)
}

/// uses a flood-fill algorithm to find all nodes connected to the given base_node.
/// the given color labels the nodes referenced in the node_colors hashmap.
fn graph_color_flood_fill(
    edge_map: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>,
    node_colors: &mut HashMap<Rc<str>, usize, FixedBuildHasher>,
    base_node: &Rc<str>,
    color: usize,
) {
    let mut node_stack: Vec<Rc<str>> = Vec::new();
    node_stack.push(base_node.clone());

    while !node_stack.is_empty() {
        let current_node: Rc<str> = node_stack.pop().unwrap();
        if !node_colors.contains_key(&current_node.clone()) {
            node_colors.insert(current_node.clone(), color);

            let node_edges = edge_map.get(&current_node).unwrap();

            // Sort edges for deterministic traversal order
            let mut sorted_in_edges: Vec<Rc<str>> = node_edges.in_edges.clone();
            let mut sorted_out_edges: Vec<Rc<str>> = node_edges.out_edges.clone();
            sorted_in_edges.sort();
            sorted_out_edges.sort();

            // Use sorted order for deterministic behavior
            node_stack.extend_from_slice(&sorted_in_edges);
            node_stack.extend_from_slice(&sorted_out_edges);
        }
    }
}

/// determines if the given graph has any cycles
fn graph_is_acyclic(node_map: &HashMap<Rc<str>, NodeEdges, FixedBuildHasher>) -> bool {
    let mut source_nodes: Vec<Rc<str>> = node_map
        .iter()
        .filter(|(_, edges)| edges.in_edges.is_empty())
        .map(|(node_name, _)| node_name.clone())
        .collect();

    // Sort source nodes for deterministic cycle detection
    source_nodes.sort();

    for source in source_nodes {
        let mut acyclic_nodes: HashSet<Rc<str>> = HashSet::new();
        let mut node_path: Vec<(Rc<str>, usize)> = Vec::new();

        node_path.push((source, 0));

        while !node_path.is_empty() {
            let (current_node, child_counter) = node_path.pop().unwrap();
            let out_edges = &node_map.get(&current_node).unwrap().out_edges;

            //if this node has more children to process and we aren't sure if it's cyclic,
            if child_counter < out_edges.len() && !acyclic_nodes.contains(&current_node) {
                //put the node back on the stack, with an updated child counter
                node_path.push((current_node, child_counter + 1));
                let child_node_name = out_edges.get(child_counter).unwrap();
                if node_path
                    .iter()
                    .any(|(node_name, _)| *node_name == *child_node_name)
                {
                    return false;
                }
                node_path.push((child_node_name.clone(), 0));
            } else {
                //if we're at the end of the node's child list, it def doesn't have any cycles.
                acyclic_nodes.insert(current_node);
            }
        }
    }

    return true;
}
