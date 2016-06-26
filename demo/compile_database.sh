#!/bin/bash
ls -d S* | awk '{print $1"\t"$1}' - > dataset_directory.txt

python ../chilin2/modules/interface/tojson.py dataset_directory.txt
