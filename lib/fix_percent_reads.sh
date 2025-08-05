#!/bin/bash
#Christie Woodside
#August 5, 2025

# Code is used to fix the percent reads aligned and percent reads unaligned values in the assemblyQC_ARGOS core table. These values were incorrect in the JSON outputs
# (now corrected). This code takes in the assemblyQC tsv and fixes the values and spits out the corrected table.

### bash fix_percent_reads.sh /Users/user/Desktop/test_assemblyQC_ARGOS.tsv /Users/user/Desktop/fixed_assemblyQC_ARGOS.tsv

input_file="$1"
output_file="${2:-fixed_output.tsv}"

if [ -z "$input_file" ]; then
  echo "Usage: $0 input_file [output_file]"
  exit 1
fi

awk -F'\t' -v OFS='\t' '
NR==1 {
  # Store header and get column indexes
  for (i=1; i<=NF; i++) {
    header[i] = $i
    colname[$i] = i
  }
  print $0
  next
}

{
  gaid = $(colname["genome_assembly_id"])
  aligned = $(colname["reads_aligned"])
  if (aligned ~ /^[0-9]+$/)
    aligned_sum[gaid] += aligned
  data[NR] = $0
  row_gaid[NR] = gaid
  row_unaligned[NR] = $(colname["reads_unaligned"])
  row_aligned[NR] = aligned
}

END {
  for (r in data) {
    gaid = row_gaid[r]
    unaligned = row_unaligned[r]
    aligned = row_aligned[r]
    total = aligned_sum[gaid] + 0

    if (unaligned ~ /^[0-9]+$/ && total > 0) {
      total_reads = unaligned + total
      percent_unaligned = sprintf("%.2f%%", (unaligned / total_reads) * 100)
      percent_aligned = sprintf("%.2f%%", (aligned / total_reads) * 100)
    } else {
      percent_unaligned = "NA"
      percent_aligned = "NA"
    }

    # Output the line, but replace percent columns
    split(data[r], fields, "\t")
    fields[colname["percent_reads_unaligned"]] = percent_unaligned
    fields[colname["percent_reads_aligned"]] = percent_aligned

    line = fields[1]
    for (i=2; i<=length(fields); i++) {
      line = line OFS fields[i]
    }
    print line
  }
}
' "$input_file" > "$output_file.unsorted"

# Sort the fixed output by genome_assembly_id (column 4)
(head -n 1 "$output_file.unsorted" && tail -n +2 "$output_file.unsorted" | sort -t $'\t' -k4,4) > "$output_file"
rm "$output_file.unsorted"

echo "✅ Fixed and sorted file written to: $output_file"
