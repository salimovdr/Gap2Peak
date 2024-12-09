#!/bin/bash
set -e

if [ -z "$1" ]; then
    echo "Enter a GSM ID as an argument."
    exit 1
fi

GSM_id="$1"
input_file="data/${GSM_id}.sam"

if [ ! -f "$input_file" ]; then
    echo "File $input_file not found."
    exit 1
fi

awk '{for(i=1;i<=NF;i++) if($i ~ /^RG:Z:/) print $i}' "$input_file" | \
    sort | uniq -c | sort -nr | head -150 | \
    awk '$2 !~ /^RG:Z:GGGGGGGG/' > "filtered_${GSM_id}.txt"

head -100 "filtered_${GSM_id}.txt" > "data/top_100_${GSM_id}_w_quantity.txt"

awk '{print $2}' "data/top_100_${GSM_id}_w_quantity.txt" > "data/top_100_${GSM_id}.txt"

rm "filtered_${GSM_id}.txt"

echo "Results saved to files top_100_${GSM_id}_w_quantity.txt and top_100_${GSM_id}.txt"
