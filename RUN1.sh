#!/bin/bash
./compile_opt.sh
./test_6 | tee out.dat
python extract_data.py
