#!/bin/bash

mkdir -p poppunk_temp

for file in $@
do
	cp ../${file} poppunk_temp/${file##*/}
done
