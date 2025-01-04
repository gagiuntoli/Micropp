#!/bin/bash

find ./src -name "*.cc" -o -name "*.h" | while read -r file; do
    echo "Formatting $file"
    clang-format-18 -i "$file"
done

