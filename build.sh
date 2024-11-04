echo "Compiling needed rust files..."
cd libs/graph_analyze/src
cargo build --release
cd ../../..
cd libs/super_source_and_sink/src
cargo build --release
cd ../../..
echo "Build Complete."
