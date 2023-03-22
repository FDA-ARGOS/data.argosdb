#! /bin/bash

SYSTEM=$(uname -s)
GAWK_CMD=$(which gawk)
REALPATH_CMD=$(which realpath)
DOCKER_CMD=$(which docker)

BOLD="\x1b[1m"
NORM="\x1b[0m"

set -e

if [ ! -f "$REALPATH_CMD" ]; then
  echo "Command dependency 'realpath' is not found/installed." >&2
  case "${SYSTEM}" in
    Darwin*)
      echo -e "To install using Homebrew (https://brew.sh/): ${BOLD}brew install coreutils${NORM}" >&2
      ;;
  esac
  exit 5
fi

if [ ! -f "$GAWK_CMD" ]; then
  echo "Command dependency 'gawk' is not found/installed." >&2
  case "${SYSTEM}" in
    Darwin*)
      echo -e "To install using Homebrew (https://brew.sh/): ${BOLD}brew install gawk${NORM}" >&2
      ;;
  esac
  exit 5
fi

if [ ! -f "$DOCKER_CMD" ]; then
  echo "Command dependency 'docker' is not found/installed." >&2
  echo -e "Visit ${BOLD}https://docs.docker.com/get-docker/${NORM} for installation." >&2
  exit 4
fi

USAGE="Usage: $0 <REF_NAME> <SIERRA_REPORT_DIR>"

ref_name="$1"
sierra_report="$2"
idx=1
for (( ; ; )); do
  if [ "$idx" -eq 1 ]; then
    test_sierra_report="$HOME/Downloads/analysis-reports"
  else
    test_sierra_report="$HOME/Downloads/analysis-reports $idx"
  fi
  idx=$((idx + 1))
  if [ ! -d "$test_sierra_report" ]; then
    break
  fi
  default_sierra_report="$test_sierra_report"
done

if [ -z "$ref_name" ]; then
  echo -n "Enter refName: "
  read ref_name
fi

test_sierra_report="$HOME/Downloads/$ref_name/analysis-reports"
if [ -d "$test_sierra_report" ]; then
  default_sierra_report="$test_sierra_report"
fi

test_sierra_report="$HOME/Downloads/$ref_name/NGS-analysis-reports"
if [ -d "$test_sierra_report" ]; then
  default_sierra_report="$test_sierra_report"
fi

while [ ! -d "$sierra_report" ]; do
  if [ -d "$default_sierra_report" ]; then
    echo -ne "Sierra-SARS2 report dir at ${BOLD}${default_sierra_report}${NORM} [Enter/type in]: "
  else
    echo -n "Sierra-SARS2 report dir [type in]: "
  fi
  read sierra_report
  if [ -z "$sierra_report" -a -d "$default_sierra_report" ]; then
    echo "  Use detected path $default_sierra_report"
    sierra_report="$default_sierra_report"
  fi
done

sierra_mutlist="$sierra_report/mutation-list"
if [ ! -d "$sierra_mutlist" ]; then
  echo "Error: missing directory $sierra_mutlist" 1>&2
  exit 6
fi

sierra_summary="$sierra_report/sequenceSummaries.csv"
if [ ! -f "$sierra_summary" ]; then
  echo "Error: missing file $sierra_summary" 1>&2
  exit 7
fi

pango_idx='$14'
is_ngs=0
if head -1 "$sierra_summary" | \grep 'Median Read Depth' > /dev/null; then
  # NGS report
  is_ngs=1
  pango_idx='$15'
fi

lower_ref_name=$(echo $ref_name | tr '[:upper:]' '[:lower:]')

isomut_file="payload/tables/isolate_mutations.d/${lower_ref_name}-isomuts.csv"

# $isomut_file should be a file
touch $isomut_file

docker run --rm -it \
	--volume=$(pwd):/covid-drdb/ \
	--volume="$(realpath "$sierra_mutlist"):/sierra-mutlist" \
	--volume="$(realpath "$isomut_file"):/output/isomuts.csv" \
 		hivdb/covid-drdb-builder:latest \
	pipenv run python -m drdb.entry extract-sierra-mutations /sierra-mutlist "/output/isomuts.csv"

if [ $is_ngs -eq 1 ]; then
  sed -i '' 's/\.codfreq[^,]*//g' "$isomut_file"
else
  sed -i '' 's/|[^,]*//g' "$isomut_file"
fi

echo -e "Create ${BOLD}${isomut_file}${NORM}"

iso_file="payload/tables/isolates.d/${lower_ref_name}-iso.csv"
echo "iso_name,var_name,site_directed,gisaid_id,genbank_accn,sra_accn,expandable" > "$iso_file"
if [ $is_ngs -eq 1 ]; then
  cat "$sierra_summary" |
    tail -n +2 |
    sort |
    gawk -vFPAT='[^,]*|"[^"]*"' "{
      name=\$1
      pango=${pango_idx}
      split(name, nameArr, \".codfreq\")
      printf(\"%s,%s,FALSE,NULL,NULL,%s,TRUE\n\", nameArr[1], pango, nameArr[1])
    }" >> "$iso_file"
else
  cat "$sierra_summary" |
    tail -n +2 |
    sort |
    gawk -vFPAT='[^,]*|"[^"]*"' "{
      name=\$1
      pango=${pango_idx}
      split(name, nameArr, \"|\")
      printf(\"%s,%s,FALSE,%s,NULL,NULL,TRUE\n\", nameArr[1], pango, nameArr[2])
    }" >> "$iso_file"
fi
# remove EPI_ISL_ initial from gisaid_id since they should be integers
sed -i '' 's/,EPI_ISL_/,/g' "$iso_file"
num_iso=$(wc -l "$iso_file" | awk '{print $1}')
num_iso=$((num_iso - 1))
echo -e "Create ${BOLD}${iso_file}${NORM} with ${BOLD}${num_iso}${NORM} isolates"

cat "$iso_file" |
  gawk -vFPAT='[^,]*|"[^"]*"' '{print $2}' |
  tail -n +2 |
  sort |
  uniq |
  while read var_name; do
    if [ -z "$(\grep -F "$var_name," payload/tables/variants.csv)" ]; then
      echo -e "  Missing ${BOLD}${var_name}${NORM}: manually add via editing payload/tables/variants.csv"
    fi
  done

for (( ; ; )); do
  breakbreak=0
  echo -n "Are these $num_iso isolates linked to one or more subjects (patients or animal models)? [Yes/No]: "
  read has_sbjiso
  case $has_sbjiso in
    [Yy]* )
      breakbreak=1

      sbjiso_file="payload/tables/subject_isolates/${lower_ref_name}-sbjiso.csv"
      echo "ref_name,subject_name,collection_date_cmp,collection_date,iso_name,iso_source,iso_culture,location,section" > "$sbjiso_file"
      if [ $is_ngs -eq 1 ]; then
        cat "$sierra_summary" |
          tail -n +2 |
          sort |
          gawk -vFPAT='[^,]*|"[^"]*"' "{
            name=\$1
            pango=${pango_idx}
            split(name, nameArr, \".codfreq\")
            printf(\"${ref_name},%s,=,NULL,%s,NP,FALSE,NULL,NULL\n\", pango, nameArr[1])
          }" >> /tmp/sbjiso.csv
      else
        cat "$sierra_summary" |
          tail -n +2 |
          sort |
          gawk -vFPAT='[^,]*|"[^"]*"' "{
            name=\$1
            pango=${pango_idx}
            split(name, nameArr, \"|\")
            printf(\"${ref_name},%s,=,%s,%s,NP,FALSE,NULL,NULL\n\", pango, nameArr[3], nameArr[1])
          }" >> /tmp/sbjiso.csv
      fi
      cat /tmp/sbjiso.csv | sort >> $sbjiso_file
      rm /tmp/sbjiso.csv
      echo -e "Create ${BOLD}${sbjiso_file}${NORM}"
      echo "  Manually edit ${sbjiso_file} to update subject info"

      sbjinf_file="payload/tables/subject_infections/${lower_ref_name}-inf.csv"
      echo "ref_name,subject_name,infection_date_cmp,infection_date,infected_var_name,location,immune_status,severity,section" > "$sbjinf_file"
      cat "$sbjiso_file" |
        tail -n +2 |
        gawk -vFPAT='[^,]*|"[^"]*"' '{print $2}' |
        uniq > /tmp/sbjs.csv
      for sbj in $(cat /tmp/sbjs.csv); do
        \grep -F "$ref_name,$sbj," $sbjiso_file |
          head -1 |
          gawk -vFPAT='[^,]*|"[^"]*"' '{printf("%s,%s,<,%s,%s,NULL,NULL,NULL,NULL\n", $1, $2, $4, $2)}' >> "$sbjinf_file"
      done
      echo -e "Create ${BOLD}${sbjinf_file}${NORM}"
      echo "  Manually edit ${sbjinf_file} to update subject info"

      break;;
    [Nn]* )
      breakbreak=1
      break;;
  esac
  if [ "$breakbreak" -eq 1 ]; then
    break
  fi
done
