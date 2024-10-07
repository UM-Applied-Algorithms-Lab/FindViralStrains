use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

// Work in progress, sorry //
fn read_seq(file: &str) -> Vec<[String; 2]> {
    let mut edges = Vec::new();
    let lines = read_lines(file).expect("Could not read lines");

    for line in lines {
        if let Ok(line) = line {
            let mut parts: Vec<_> = line.split_whitespace().collect();
            parts.remove(0); // Remove the name
            let (mut nlist, mut slist): (Vec<String>, Vec<String>) = parts
                .iter()
                .map(|&j| {
                    let len = j.len();
                    (j[..len - 1].to_string(), j[len - 1..].to_string())
                })
                .unzip();

            for k in 0..nlist.len() - 1 {
                let edge = if k >= 0 {
                    [nlist[k].clone() + &slist[k], nlist[k + 1].clone() + &slist[k + 1]]
                } else {
                    [
                        nlist[k + 1].clone() + if slist[k + 1] == "-" { "+" } else { "-" },
                        nlist[k].clone() + if slist[k] == "-" { "+" } else { "-" },
                    ]
                };
                edges.push(edge);
            }
        }
    }
    edges
}

fn revcomp(seq: &str) -> String {
    let mut rcseq = seq
        .replace("A", "t")
        .replace("C", "g")
        .replace("T", "a")
        .replace("G", "c")
        .to_uppercase();
    rcseq = rcseq.chars().rev().collect::<String>();
    rcseq
}

fn read_seg(file: &str, edges: Vec<[String; 2]>, k: usize) -> HashMap<String, String> {
    let mut node_seg = HashMap::new();
    if let Ok(lines) = read_lines(file) {
        for line in lines {
            if let Ok(line) = line {
                let parts: Vec<_> = line.split_whitespace().collect();
                node_seg.insert(parts[0].to_string() + "+", parts[1].to_string());
                node_seg.insert(parts[0].to_string() + "-", revcomp(parts[1]));
            }
        }
    }

    let mut edge_mers = HashMap::new();
    for edge in edges {
        let mer = format!(
            "{}{}",
            &node_seg[&edge[0]][node_seg[&edge[0]].len() - k..],
            &node_seg[&edge[1]][..k]
        );
        edge_mers.insert(edge.join(" "), mer);
    }

    edge_mers
}

fn write_graph(file: &str, edges: Vec<[String; 2]>, edge_mers: HashMap<String, String>) {
    let mut nodes = HashSet::new();
    for e in edges.iter() {
        nodes.insert(&e[0]);
        nodes.insert(&e[1]);
    }

    let mut output = File::create(file).expect("Could not create output file");
    writeln!(output, "#graph 1").expect("Failed to write");
    writeln!(output, "{}", nodes.len()).expect("Failed to write");

    for (e, m) in edge_mers.iter() {
        writeln!(output, "{} {}", e, m).expect("Failed to write");
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} -k <kvalue> -c <cuttlefishprefix> -o <output>", args[0]);
        return;
    }

    let kvalue: usize = args[2].parse().expect("Invalid K value");
    let cuttlefish_prefix = &args[4];
    let output_file = &args[6];

    println!("reading seq");
    let edges = read_seq(&(cuttlefish_prefix.clone() + ".cf_seq"));
    println!("reading seg");
    let edge_mers = read_seg(&(cuttlefish_prefix.clone() + ".cf_seg"), edges, kvalue);
    println!("writing graph");
    write_graph(output_file, edges, edge_mers);
}
