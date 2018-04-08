#!/bin/bash
./compile_opt.sh
./test_6 | tee out.dat
./plot_output.sh out.dat
python extract_data.py y
