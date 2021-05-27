#!/bin/bash
set -e
(cd src_basic && make clean && make)
(cd src_comb && make clean && make)
(cd src_graph && make clean && make)
(cd src_matrix && make clean && make)
(cd src_number && make clean && make)
(cd src_sat && make clean && make)
