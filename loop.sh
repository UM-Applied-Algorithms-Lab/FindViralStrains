
#!/bin/bash

for sample_dir in test/data/NoRefTest/dbg/*/; do
    sample=$(basename "$sample_dir")
    input_file="test/data/NoRefTest/dbg/$sample/out.dbg"
    
    for i in {0..1}; do
        # Paths for pruning
        pruned_dir="test/data/NoRefTest/dbg/$sample/pruned_$i"
	pruned_file="$pruned_dir/out.dbg"
	pruned_ref_file="$pruned_dir/ref_out.dbg"
   	mkdir -p "$pruned_dir"

        # Run pruning
        echo "[$sample] Pruning with num_prune=$i"
	python3 libs/prune/add_reference.py "$pruned_file" "reference_genomes/covid19ref.fasta" "$pruned_ref_file"

        # Paths for subgraph analysis
        subgraph_dir="$pruned_dir/out.dbg_subgraphs"
        stats_dir="dbg/$sample"
        stats_file="$stats_dir/graph_stats_$i.txt"
        mkdir -p "$subgraph_dir" "$stats_dir"

        # Run graph_analyzer
        echo "[$sample] Analyzing pruned graph (num_prune=$i)"
        target/release/graph_analyzer \
		 --dbg-file-name "$pruned_ref_file" --node-percent-cutoff 0 > "$stats_file"

    done
done

