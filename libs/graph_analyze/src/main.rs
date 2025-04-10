use clap::Parser;
use colored::Colorize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Write};
use std::path::Path;
use std::rc::Rc;

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
}

//struct to hold the analysis data for a graph or subgraph
struct GraphAnalysisData {
    num_nodes: usize,
    num_edges: usize,
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
            "{}\t\t{}\t\t{}\t{}\t\t{}\t\t{}\n{}\t\t{}\t\t{}\t\t{}\t\t{}\t\t{}",
            "nodes".blue(),
            "edges".blue(),
            "subgraphs".blue(),
            "sources".blue(),
            "sinks".blue(),
            "acyclic?".blue(),
            self.num_nodes,
            self.num_edges,
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
    );
}

/// given a list of subgraphs, writes them to files for further use in the pipeline.
/// writes subgraph source, sink, and .mg files to a subdirectory named after the input .mg file src
///
/// significant_subgraph_list should be the list of all subgraphs to print
fn write_subgraph_files(
    significant_subgraph_list: Vec<&HashMap<Rc<str>, NodeEdges>>,
    base_file_name: &String,
    main_graph_label: &String,
    edge_info: &HashMap<(Rc<str>, Rc<str>), (Rc<str>, Rc<str>)>,
    output_dir: Option<String>,
) {
    for (subgraph_idx, subgraph) in significant_subgraph_list.iter().enumerate() {
        let subgraph_sub_dir = base_file_name.to_string() + "_subgraphs";
        let subgraph_directory_name = match output_dir {
            Some(ref dir) => Path::new(dir).join(Path::new(&subgraph_sub_dir)),
            None => Path::new(&subgraph_sub_dir).to_path_buf(),
        };

        std::fs::create_dir_all(&subgraph_directory_name).expect(
            "could not create subgraph directory"
        );
        let subgraph_idx_string = subgraph_idx.to_string();
        // let subgraph_mg_file_src = base_file_name.to_owned() + "." + &subgraph_idx_string;
        let (sources, sinks) = make_source_sink_lists(&subgraph);

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

        for (from_node, edges) in *subgraph {
            for to_node in &edges.out_edges {
                match edge_info.get(&(from_node.clone(), to_node.clone())) {
                    Some(result) => {
                        let (count, kmer) = result;
                        subgraph_mg_file
                            .write_fmt(format_args!(
                                "{}\t{}\t{}\t{}\n",
                                from_node,
                                to_node,
                                count,
                                kmer
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
    main_graph: &HashMap<Rc<str>, NodeEdges>,
    subgraph_list: &Vec<HashMap<Rc<str>, NodeEdges>>,
    subgraph_display_type: SubgraphDisplayType,
    stats_output_file: &Option<String>,
) -> std::io::Result<()> {
    let mut output = String::new();
    
    output.push_str(&format!(
        "Main Graph Stats:\n{}\n",
        make_graph_stats(main_graph, subgraph_list.len())
    ));

    match subgraph_display_type {
        SubgraphDisplayType::All => {
            for (subgraph_idx, subgraph) in subgraph_list.iter().enumerate() {
                output.push_str(&format!(
                    "Subgraph {}:\n{}\n",
                    subgraph_idx,
                    make_graph_stats(subgraph, 0)
                ));
            }
        }

        SubgraphDisplayType::OnlySignificant(cutoff) => {
            let significant_subgraph_list = make_significant_subgraph_list(&subgraph_list, cutoff);
            for (subgraph_idx, subgraph) in significant_subgraph_list.iter().enumerate() {
                output.push_str(&format!(
                    "Subgraph {}:\n{}\n",
                    subgraph_idx,
                    make_graph_stats(subgraph, 0)
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
    graph: &HashMap<Rc<str>, NodeEdges>,
    num_subgraphs: usize,
) -> GraphAnalysisData {
    let (sources, sinks) = make_source_sink_lists(&graph);
    GraphAnalysisData {
        num_nodes: graph.len(),
        num_edges: graph.iter().map(|(_, edges)| edges.out_edges.len()).sum(),
        num_disconnected_subgraphs: num_subgraphs,
        sources,
        sinks,
        is_acyclic: graph_is_acyclic(&graph),
    }
}

/// generates the list of subgraphs. This is a new graph list, so it increaes memory, but avoids memory issues and the
/// borrow checker (and lifetimes)
fn make_subgraph_list(
    main_graph: &HashMap<Rc<str>, NodeEdges>,
) -> Vec<HashMap<Rc<str>, NodeEdges>> {
    let mut node_colors: HashMap<Rc<str>, usize> = HashMap::new();
    let mut subgraph_list: Vec<HashMap<Rc<str>, NodeEdges>> = Vec::new();

    // flood fills the graph to find all connected nodes
    for (node_idx, (node_name, _)) in main_graph.iter().enumerate() {
        graph_color_flood_fill(main_graph, &mut node_colors, node_name, node_idx);
    }

    // finds the number of colors that were used to color all nodes of the graph, and therefore,
    // finds the set (and count) of subgraphs
    let node_color_set: HashSet<usize> = node_colors
        .iter()
        .map(|(_, node_color)| node_color.clone())
        .collect();

    // generates separate hashmaps for each subgraph
    for node_color in node_color_set {
        let subgraph: HashMap<Rc<str>, NodeEdges> = node_colors
            .iter()
            .filter(|(_, color)| node_color == **color)
            .map(|(node_name, _)| main_graph.get_key_value(node_name).unwrap())
            .map(|(node_name, node_edges)| (node_name.clone(), node_edges.clone()))
            .collect();
        subgraph_list.push(subgraph);
    }
    return subgraph_list;
}

/// using a percent cutoff, filters the list of subgraphs to only those who have a large enough share of nodes.
/// "large enough share" is defined as at least 'node_percent_cutoff' percent of the total nodes of the full graph.
fn make_significant_subgraph_list(
    subgraph_list: &Vec<HashMap<Rc<str>, NodeEdges>>,
    node_percent_cutoff: f32,
) -> Vec<&HashMap<Rc<str>, NodeEdges>> {
    let num_total_nodes: usize = subgraph_list.iter().map(|subgraph| subgraph.len()).sum();
    let min_nodes_for_significant_subgraph: usize =
        ((num_total_nodes as f32) * node_percent_cutoff) as usize;
    return subgraph_list
        .iter()
        .filter(|subgraph| subgraph.len() >= min_nodes_for_significant_subgraph)
        .collect();
}

/// reads in the mer-graph from the file_path, and generates the full de bruijn graph as a hashmap
/// of outgoing edges and a hashmap of the kmers labeling the edges
fn make_main_graph(
    file_path: &Path,
) -> Result<(
    HashMap<Rc<str>, NodeEdges>,
    HashMap<(Rc<str>, Rc<str>), (Rc<str>, Rc<str>)>,
    String,
)> {
    let file = File::open(file_path)?;
    let file_reader = BufReader::new(file);

    let mut main_graph: HashMap<Rc<str>, NodeEdges> = HashMap::new();
    let mut edge_info: HashMap<(Rc<str>, Rc<str>), (Rc<str>, Rc<str>)> = HashMap::new();
    let mut lines = file_reader.lines();
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

        let count: Rc<str> = Rc::from(split_line.next().expect( // TODO change to int
            "encountered incorrectly formatted line in input file: \
            could not parse second node name",
        ));

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
fn make_source_sink_lists(edge_map: &HashMap<Rc<str>, NodeEdges>) -> (Vec<Rc<str>>, Vec<Rc<str>>) {
    // let end_nodes: HashSet<&String> = out_edge_map.values().flatten().collect();

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
    edge_map: &HashMap<Rc<str>, NodeEdges>,
    node_colors: &mut HashMap<Rc<str>, usize>,
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
            node_stack.extend_from_slice(&node_edges.in_edges);
            node_stack.extend_from_slice(&node_edges.out_edges);
        }
    }
}

/// determines if the given graph has any cycles
fn graph_is_acyclic(node_map: &HashMap<Rc<str>, NodeEdges>) -> bool {
    let source_nodes: Vec<Rc<str>> = node_map
        .iter()
        .filter(|(_, edges)| edges.in_edges.is_empty())
        .map(|(node_name, _)| node_name.clone())
        .collect();

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