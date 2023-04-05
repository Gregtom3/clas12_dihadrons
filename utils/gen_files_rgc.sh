#!/bin/bash

# Set the search directory
search_dir="/volatile/clas12/rg-c/production"

# Set the output file
output_file="sidisdvcs_files.txt"

# Initialize the output file
echo -e "Basename\tNumber\tType\tDate of Creation\tPath" > "${output_file}"

# Find files with the basename "sidisdvcs_#.hipo" and sort by path and creation date
find "${search_dir}" -type f -name "sidisdvcs_*.hipo" -printf "%C@\t%p\n" | sort -k2 | cut -f 2 | while read -r file; do
  # Extract the basename
  basename=$(basename "${file}")

  # Extract the number from the basename
  number=$(echo "${basename}" | sed -n -E 's/sidisdvcs_([0-9]+).hipo/\1/p')

  # Get the creation date of the file
  creation_date=$(date -u -d @$(stat -c %Y "${file}") +"%m/%d/%Y")

  # Get the full path to the file
  path=$(echo "${file}")

  # Determine the type (TBT or HBT) based on the directory name
  type=$(echo "${file}" | grep -oE 'TBT|HBT')

  # Write the output to the file
  echo -e "${basename}\t${number}\t${type}\t${creation_date}\t${path}" >> "${output_file}"

  # Check if the next file is in the same directory, if not, add an empty row
  next_file=$(find "${search_dir}" -type f -name "sidisdvcs_*.hipo" -printf "%C@\t%p\n" | sort -k2 | cut -f 2 | grep -A1 "${file}" | tail -1)
  current_dir=$(dirname "${file}")
  next_dir=$(dirname "${next_file}")

  if [[ "${current_dir}" != "${next_dir}" ]]; then
    echo "" >> "${output_file}"
  fi
done
