#! /bin/bash

folder=$1

if [[ "$1" == *"SOMD"* ]]; then
echo "making sure there is a header as it's SOMD files..."
# search through the folder for all perts
for trans_dir in $(find $1 -name '*~*'); do
python add_header_simfile.py $trans_dir
done
fi

echo "extracting output for $folder"
python ../python/scripts/extract_output.py $folder

echo "done."

