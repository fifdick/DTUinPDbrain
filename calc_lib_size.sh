#!/usr/bin/env bash

for filename in ./rawData/*/quant.sf; do
    # Will print */ if no directories are available
    sample=$(echo "$filename" | cut -d'/' -f3)
#    echo "$sample"
    awk -v s="$sample" '{sum = sum + int($5)} END {print s "\t" sum}' "$filename"
done
