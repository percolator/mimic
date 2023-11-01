#!/bin/bash

src_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"/src
build_dir="$( cd && pwd )"/build/mimic
#build_dir="$(mktemp -d --tmpdir mimic_build_XXXX)";

mkdir -p $build_dir
cd $build_dir
cmake -DCMAKE_INSTALL_PREFIX=/usr $src_dir && make && sudo make install
