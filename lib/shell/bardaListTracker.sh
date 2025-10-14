#!/usr/bin/env bash
# BARDA organism finder — writes to .txt file, lists unmatched names,
# and reports total unique assemblies across all BARDA organisms.
set -euo pipefail

#change input parameters here
BARDA_LIST="${1:-bardaList.txt}"
DATA_FILE="${2:-/Users/christiewoodside/Desktop/biosampleMeta_ARGOS_extended.tsv}" #update manually here
#DATA_FILE="${2:-/Users/christiewoodside/Desktop/biosampleMeta_ARGOS.tsv}"

#--------------------------------------------------------------------------------

ASM_COL="${ASM_COL:-11}"          # assembly accession column
STRIP_VERSION="${STRIP_VERSION:-0}"  # strip .version suffix?
OUT_FILE="${OUT_FILE:-barda_counts_$(date +%Y%m%d_%H%M%S).txt}"

TMP_ASM_LIST=$(mktemp)  # collects all assemblies seen
: > "$OUT_FILE"

while IFS= read -r line; do
  [[ -z "${line// }" ]] && continue

  {
    printf "Total number of entries for: %s:\n" "$line"

    # === TOTAL COUNT ===
    awk -F '\t' -v q="$line" '
      function norm(s, t){ t=tolower(s); gsub(/[[:space:]]+/, " ", t); sub(/^ +/,"",t); sub(/ +$/,"",t); return t }
      function denoise_taxon(s){ sub(/ virus$/, "", s); sub(/ bacterium$/, "", s); return s }
      function strip_parens(s, p,q,left,right){
        while ((p=index(s,"("))>0){ left=substr(s,1,p-1); q=index(substr(s,p+1),")"); if(!q)break; s=left substr(s,p+q+1) }
        return s
      }
      function inner_parens(s,p,q){ p=index(s,"("); if(!p)return""; q=index(substr(s,p+1),")"); if(!q)return""; return substr(s,p+1,q-1) }
      function mk_aliases(src,raw,out,inner){ raw=src; out=strip_parens(raw); inner=inner_parens(raw);
        aliases[0]=denoise_taxon(norm(out)); aliases[1]=denoise_taxon(norm(inner)); }
      BEGIN{ mk_aliases(q) }
      { org=denoise_taxon(norm($1));
        if ((aliases[0]!=""&&(org==aliases[0]||index(org,aliases[0])>0))||
            (aliases[1]!=""&&(org==aliases[1]||index(org,aliases[1])>0))) c++ }
      END{ print c+0 }
    ' "$DATA_FILE"

    # === UNIQUE ASSEMBLIES ===
    printf "\tUnique assemblies:\n"
    awk -F '\t' -v q="$line" -v ASM_COL="$ASM_COL" -v STRIP_VERSION="$STRIP_VERSION" -v tmp="$TMP_ASM_LIST" '
      function norm(s,t){ t=tolower(s); gsub(/[[:space:]]+/, " ", t); sub(/^ +/,"",t); sub(/ +$/,"",t); return t }
      function denoise_taxon(s){ sub(/ virus$/, "", s); sub(/ bacterium$/, "", s); return s }
      function strip_parens(s,p,q,left,right){
        while ((p=index(s,"("))>0){ left=substr(s,1,p-1); q=index(substr(s,p+1),")"); if(!q)break; s=left substr(s,p+q+1) }
        return s
      }
      function inner_parens(s,p,q){ p=index(s,"("); if(!p)return""; q=index(substr(s,p+1),")"); if(!q)return""; return substr(s,p+1,q-1) }
      function mk_aliases(src,raw,out,inner){ raw=src; out=strip_parens(raw); inner=inner_parens(raw);
        aliases[0]=denoise_taxon(norm(out)); aliases[1]=denoise_taxon(norm(inner)); }
      function asm_key(a,k){ k=a;
        if (k ~ /^GC[AF]_[0-9]+(\.[0-9]+)?$/ && STRIP_VERSION==1) sub(/\.[0-9]+$/, "", k);
        return k; }
      BEGIN{ mk_aliases(q) }
      {
        org=denoise_taxon(norm($1))
        if ((aliases[0]!=""&&(org==aliases[0]||index(org,aliases[0])>0))||
            (aliases[1]!=""&&(org==aliases[1]||index(org,aliases[1])>0))){
          if (ASM_COL<=NF){
            k=asm_key($ASM_COL)
            if (k!="" && !(k in seen)){ seen[k]=1; u++; print k >> tmp }
          }
        }
      }
      END{ print u+0 }
    ' "$DATA_FILE"

    printf "\n"
  } >> "$OUT_FILE"

done < "$BARDA_LIST"

# === UNMATCHED ORGANISMS ===
{
  printf "=== Not matched in BARDA list (unique organism names from biosample table) ===\n"
  awk -F '\t' -v BL="$BARDA_LIST" '
    function norm(s,t){ t=tolower(s); gsub(/[[:space:]]+/, " ", t); sub(/^ +/,"",t); sub(/ +$/,"",t); return t }
    function denoise_taxon(s){ sub(/ virus$/, "", s); sub(/ bacterium$/, "", s); return s }
    function strip_parens(s,p,q,left,right){
      while ((p=index(s,"("))>0){ left=substr(s,1,p-1); q=index(substr(s,p+1),")"); if(!q)break; s=left substr(s,p+q+1) }
      return s
    }
    function inner_parens(s,p,q){ p=index(s,"("); if(!p)return""; q=index(substr(s,p+1),")"); if(!q)return""; return substr(s,p+1,q-1) }
    BEGIN{
      while((getline l<BL)>0){ if(l~/^[[:space:]]*$/)continue;
        out=denoise_taxon(norm(strip_parens(l))); inr=denoise_taxon(norm(inner_parens(l)));
        if(out!=""&&!(out in seenA)){ seenA[out]=1; aliases[++n]=out }
        if(inr!=""&&!(inr in seenA)){ seenA[inr]=1; aliases[++n]=inr }
      }
      close(BL)
    }
    {
      if(NR==1 && tolower($1)~/organism|species|taxon/) next
      raw=$1; org=denoise_taxon(norm(raw)); matched=0
      for(i=1;i<=n;i++){ a=aliases[i]; if(a!=""&&(org==a||index(org,a)>0)){ matched=1; break } }
      if(!matched && !(raw in outed)){ outed[raw]=1; print raw }
    }
  ' "$DATA_FILE" | sort -fu
} >> "$OUT_FILE"

# === GRAND TOTAL OF UNIQUE ASSEMBLIES ===
{
  printf "\n=== Total unique assemblies mapping to BARDA list ===\n"
  sort -u "$TMP_ASM_LIST" | wc -l | awk '{print $1}'
} >> "$OUT_FILE"

rm -f "$TMP_ASM_LIST"
printf "Wrote results to: %s\n" "$OUT_FILE"




# ----------------------------------------------------
#og stuff
# #!/bin/bash
# # Script to look through our processed data to find those requested by BARDA. This file looks through the assemblyQC_ARGOS" file.
# # Christie's updates on October 14

# while read line; do
#         echo "Total number of entries for: $line:"
#         #cat /data/shared/argosdb/generated/datasets/reviewed/assemblyQC_ARGOS_extended.tsv | awk -F '\t' -v name="$line" '$1 ~ name' | wc -l

#         cat /Users/christiewoodside/Desktop/biosampleMeta_ARGOS_extended.tsv | awk -F '\t' -v name="$line" '$1 ~ name' | wc -l
#         echo -e ' \t '"Unique assembles:"
#         #echo -ne ' \t ' && cat /data/shared/argosdb/generated/datasets/reviewed/assemblyQC_ARGOS_extended.tsv |  awk -F '\t' -v name="$line" '$1 ~ name' | awk -F '\t' '{ print $4 }' | sort | uniq | wc -l

#         echo -ne ' \t ' && cat /Users/christiewoodside/Desktop/biosampleMeta_ARGOS_extended.tsv |  awk -F '\t' -v name="$line" '$1 ~ name' | awk -F '\t' '{ print $4 }' | sort | uniq | wc -l
# done < bardaList.txt