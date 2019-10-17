#!/bin/bash

mkdir poppunk_export

for dir in $@
do
	cp ${dir}/K?_core_NJ.nwk ${dir}/K?_microreact_clusters.csv ${dir}/K?_perplexity*_accessory_tsne.dot poppunk_export
done

zip poppunk_export.zip poppunk_export/*
