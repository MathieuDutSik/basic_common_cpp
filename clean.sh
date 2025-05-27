#!/bin/bash
set -e
(cd src_basic && make clean)
(cd src_comb && make clean)
(cd src_graph && make clean)
(cd src_matrix && make clean)
(cd src_number && make clean)
(cd src_sat && make clean)
