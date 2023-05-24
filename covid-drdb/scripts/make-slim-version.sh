#!/bin/bash

DB_FILE=$1

sqlite3 $DB_FILE <<EOF
DELETE FROM subject_plasma;
DELETE FROM rx_potency;
DELETE FROM rx_fold;
DELETE FROM treatments;
DELETE FROM ref_invivo;
DELETE FROM ref_isolate_pairs;
DELETE FROM ref_unpaired_isolates;
DELETE FROM antibody_articles;
DELETE FROM dms_escape_results AS d WHERE NOT EXISTS (
  SELECT 1 FROM amino_acid_prevalence p
  WHERE
    proportion > 0.00001 AND
    d.gene = p.gene AND
    d.position = p.position AND
    d.amino_acid = p.amino_acid
) OR escape_score < 0.1;
DELETE FROM dms_ace2_binding AS d WHERE NOT EXISTS (
  SELECT 1 FROM dms_escape_results p
  WHERE
    d.gene = p.gene AND
    d.position = p.position AND
    d.amino_acid = p.amino_acid
);
DELETE FROM amino_acid_prevalence;
DELETE FROM isolate_mutations WHERE gene != 'S';
DELETE FROM invivo_selection_results WHERE gene != 'S';

DELETE FROM articles AS R
  WHERE
    NOT EXISTS (
      SELECT 1 FROM susc_summary S
      WHERE R.ref_name = S.ref_name
    ) AND NOT EXISTS (
      SELECT 1 FROM variant_status VS
      WHERE R.ref_name = VS.ref_name
    ) AND NOT EXISTS (
      SELECT 1 FROM invitro_selection_results IM
      WHERE R.ref_name = IM.ref_name
    ) AND NOT EXISTS (
      SELECT 1 FROM invivo_selection_results IM
      WHERE R.ref_name = IM.ref_name
    ) AND NOT EXISTS (
      SELECT 1 FROM dms_escape_results DM
      WHERE R.ref_name = DM.ref_name
    );

DELETE FROM antibody_targets AS A
  WHERE NOT EXISTS (
    SELECT 1 FROM susc_summary S
    WHERE
      S.aggregate_by = 'antibody:indiv' AND
      A.ab_name = S.antibody_names
  ) AND NOT EXISTS (
    SELECT 1 FROM rx_antibodies RXMAB, invitro_selection_results IM
    WHERE
      IM.ref_name = RXMAB.ref_name AND
      IM.rx_name = RXMAB.rx_name AND
      RXMAB.ab_name = A.ab_name
  ) AND NOT EXISTS (
    SELECT 1 FROM rx_antibodies RXMAB, dms_escape_results DM
    WHERE
      DM.ref_name = RXMAB.ref_name AND
      DM.rx_name = RXMAB.rx_name AND
      RXMAB.ab_name = A.ab_name
  );

DELETE FROM antibodies AS A
  WHERE NOT EXISTS (
    SELECT 1 FROM susc_summary S
    WHERE
      S.aggregate_by = 'antibody:indiv' AND
      A.ab_name = S.antibody_names
  ) AND NOT EXISTS (
    SELECT 1 FROM rx_antibodies RXMAB, invitro_selection_results IM
    WHERE
      IM.ref_name = RXMAB.ref_name AND
      IM.rx_name = RXMAB.rx_name AND
      RXMAB.ab_name = A.ab_name
  ) AND NOT EXISTS (
    SELECT 1 FROM rx_antibodies RXMAB, dms_escape_results DM
    WHERE
      DM.ref_name = RXMAB.ref_name AND
      DM.rx_name = RXMAB.rx_name AND
      RXMAB.ab_name = A.ab_name
  );

DELETE FROM susc_summary;
DELETE FROM resistance_mutation_attributes;
VACUUM
EOF
