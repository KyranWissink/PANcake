#!/bin/bash

for file in *.fas; do
    output="processed/${file%.fas}.fa"
    awk '{ if (substr($0, 1, 1) == ">") { if (seq) print seq; print $0; seq=""; } else { seq = seq $0 } } END { if (seq) print seq }' "$file" | fold -w 60 > "$output"
    echo "Processed $file -> $output"
done
