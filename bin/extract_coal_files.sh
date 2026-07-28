#!/bin/bash

#Define files and paths
COAL_FILE=$1
#COAL_FILE=/homes/users/xmarti/scratch/gcat/relate/ALL_RESULTS/1000GP_Phase3_mask_prene.coal
POP=$2
#POP=FIN
BASE_DIR=$3
OUTDIR=${BASE_DIR}/${POP}/coal
#OUTDIR=/homes/users/xmarti/scratch/gcat/relate/${POP}/coal

#Make necessary directory
mkdir -p $OUTDIR

# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <POP>"
    exit 1
fi

# Read the input file
mapfile -t lines < $COAL_FILE

# Extract the header POPs
headers=(${lines[0]})

# Find the position of the specified POP in the headers
POP_position=-1
for i in "${!headers[@]}"; do
    if [ "${headers[i]}" == "$POP" ]; then
        POP_position=$i
        break
    fi
done

# Check if the provided POP is in the header
if [ "$POP_position" -eq -1 ]; then
    echo "Element $POP not found in the input file."
    exit 1
fi

# Output filename
output_filename="${OUTDIR}/${POP}.coal"

# Write the first row as 0
echo "0" > "$output_filename"

# Write the second row from the original file
echo "${lines[1]}" >> "$output_filename"

# Find and write the corresponding row for the specified POP
for line in "${lines[@]:2}"; do
    first_column=$(echo "$line" | awk '{print $1}')
    second_column=$(echo "$line" | awk '{print $2}')

    if [ "$first_column" -eq "$POP_position" ] && [ "$first_column" -eq "$second_column" ]; then
        echo "$line" >> "$output_filename"
        break
    fi
done

# Overwrite the first two columns of the third row with "0 0"
sed -i '3s/^[0-9]* [0-9]*/0 0/' "$output_filename"

echo "File '$output_filename' generated for population '$POP'." > ${POP}_coal.report

