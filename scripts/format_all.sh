#!/bin/bash

find ./src ./include -name "*.cpp" -o -name "*.hpp" -o -name "*.h" | while read -r file; do
    echo "Formatting $file"
    clang-format-18 -i "$file"
done

