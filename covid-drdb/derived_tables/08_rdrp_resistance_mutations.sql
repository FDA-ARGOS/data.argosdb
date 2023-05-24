SELECT drug_name, abbreviation_name
  INTO TABLE approved_drugs
  FROM compounds
  WHERE abbreviation_name IN ('RDV', 'MOL');

INSERT INTO resistance_mutation_attributes
SELECT
  im.gene,
  position,
  amino_acid,
  'FOLD:' || drug.drug_name AS col_name,
  ROUND(CAST(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY fold) AS NUMERIC), 1) AS col_value
  FROM isolate_mutations im
    JOIN susc_results s ON
      im.iso_name = s.iso_name
    JOIN isolates iso ON
      im.iso_name = iso.iso_name
    JOIN isolate_pairs ip ON
      ip.gene = 'RdRP' AND
      s.control_iso_name = ip.control_iso_name AND
      s.iso_name = ip.iso_name
    JOIN rx_compounds rxdrug ON
      s.ref_name=rxdrug.ref_name AND
      s.rx_name=rxdrug.rx_name
    JOIN approved_drugs drug ON
      rxdrug.drug_name = drug.drug_name
    WHERE
      im.gene = 'RdRP' AND
      ip.num_mutations = 1 AND
      (SELECT COUNT(*)
        FROM rx_compounds rxdrug2
        WHERE
          s.ref_name=rxdrug2.ref_name AND
          s.rx_name=rxdrug2.rx_name
      ) = 1
  GROUP BY im.gene, position, amino_acid, col_name
  ORDER BY im.gene, position, amino_acid;

SELECT
  im.gene, position, amino_acid, s.ref_name
  INTO _drm_articles
  FROM isolate_mutations im
    JOIN susc_results s ON
      im.iso_name = s.iso_name
    JOIN isolates iso ON
      im.iso_name = iso.iso_name
    JOIN isolate_pairs ip ON
      ip.gene = 'RdRP' AND
      s.control_iso_name = ip.control_iso_name AND
      s.iso_name = ip.iso_name
    JOIN rx_compounds rxdrug ON
      s.ref_name=rxdrug.ref_name AND
      s.rx_name=rxdrug.rx_name
    JOIN approved_drugs drug ON
      rxdrug.drug_name = drug.drug_name
    WHERE
      im.gene = 'RdRP' AND
      ip.num_mutations = 1 AND
      (SELECT COUNT(*)
        FROM rx_compounds rxdrug2
        WHERE
          s.ref_name=rxdrug2.ref_name AND
          s.rx_name=rxdrug2.rx_name
      ) = 1
  GROUP BY im.gene, position, amino_acid, s.ref_name;

INSERT INTO resistance_mutation_attributes
SELECT
  gene, position, amino_acid,
  'INVIVO' AS col_name,
  SUM(count) AS col_value
  FROM invivo_selection_results sel
  WHERE
    gene = 'RdRP' AND
    amino_acid != 'stop' AND
    EXISTS (
      SELECT 1 FROM subject_treatments sbjrx
        JOIN rx_compounds rxdrug ON
          sbjrx.ref_name = rxdrug.ref_name AND
          sbjrx.rx_name = rxdrug.rx_name
        JOIN approved_drugs drug ON
          rxdrug.drug_name = drug.drug_name
      WHERE
        sbjrx.ref_name = sel.ref_name AND
        sbjrx.subject_name = sel.subject_name AND
        sbjrx.start_date < sel.appearance_date
    )
  GROUP BY gene, position, amino_acid
  ORDER BY gene, position, amino_acid;

INSERT INTO _drm_articles
SELECT
  gene, position, amino_acid, ref_name
  FROM invivo_selection_results sel
  WHERE
    gene = 'RdRP' AND
    amino_acid != 'stop' AND
    EXISTS (
      SELECT 1 FROM subject_treatments sbjrx
        JOIN rx_compounds rxdrug ON
          sbjrx.ref_name = rxdrug.ref_name AND
          sbjrx.rx_name = rxdrug.rx_name
        JOIN approved_drugs drug ON
          rxdrug.drug_name = drug.drug_name
      WHERE
        sbjrx.ref_name = sel.ref_name AND
        sbjrx.subject_name = sel.subject_name AND
        sbjrx.start_date < sel.appearance_date
    )
  GROUP BY gene, position, amino_acid, ref_name
ON CONFLICT DO NOTHING;

INSERT INTO resistance_mutation_attributes
SELECT
  gene, position, amino_acid,
  'INVITRO' AS col_name,
  COUNT(*) AS col_value
  FROM invitro_selection_results
  WHERE gene = 'RdRP' AND amino_acid != 'stop'
  GROUP BY gene, position, amino_acid
  ORDER BY gene, position, amino_acid;

INSERT INTO _drm_articles
SELECT
  gene, position, amino_acid, ref_name
  FROM invitro_selection_results
  WHERE gene = 'RdRP' AND amino_acid != 'stop'
  GROUP BY gene, position, amino_acid, ref_name
ON CONFLICT DO NOTHING;

INSERT INTO resistance_mutation_attributes
SELECT
  gene, position, amino_acid,
  'PREVALENCE' AS col_name,
  proportion AS col_value
  FROM amino_acid_prevalence
  WHERE gene = 'RdRP' AND ref_name = 'Martin21'
  ORDER BY gene, position, amino_acid;

INSERT INTO _drm_articles
SELECT
  gene, position, amino_acid, ref_name
  FROM amino_acid_prevalence
  WHERE gene = 'RdRP' AND ref_name = 'Martin21'
  ORDER BY gene, position, amino_acid, ref_name
ON CONFLICT DO NOTHING;

INSERT INTO resistance_mutations
SELECT
  gene, position, amino_acid
FROM resistance_mutation_attributes rma
  WHERE
    gene = 'RdRP' AND (
      (col_name LIKE 'FOLD:%' AND
       col_value::DECIMAL >= 2.5) OR
      (col_name = 'INVITRO' AND
       col_value::DECIMAL > 1) OR
      (col_name = 'INVIVO' AND
       col_value::DECIMAL > 1)
    ) AND
    NOT EXISTS (
      SELECT 1
      FROM ignore_mutations igm
      WHERE
        igm.gene = rma.gene AND
        igm.position = rma.position AND
        igm.amino_acid = rma.amino_acid
    )
  GROUP BY gene, position, amino_acid;

INSERT INTO resistance_mutation_attributes
SELECT
  vc.gene, vc.position, vc.amino_acid,
  'VARCONS' AS col_name,
  '1' AS col_value
  FROM variant_consensus vc
    JOIN resistance_mutations rm ON
      vc.gene = rm.gene AND
      vc.position = rm.position AND
      vc.amino_acid = rm.amino_acid
  WHERE
    vc.gene = 'RdRP' AND
    vc.var_name IN (
      'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon',
      'Zeta', 'Eta', 'Theta', 'Iota', 'Kappa',
      'Lambda', 'Mu', 'Omicron/BA.1', 'Omicron/BA.2', 'Omicron/BA.3',
      'Omicron/BA.4', 'Omicron/BA.5'
    )
  GROUP BY vc.gene, vc.position, vc.amino_acid, col_name, col_value
  ORDER BY vc.gene, vc.position, vc.amino_acid, col_name, col_value;

DELETE FROM resistance_mutation_attributes rma
  WHERE
    gene = 'RdRP' AND
    NOT EXISTS (
      SELECT 1 FROM resistance_mutations rm
      WHERE
        rm.gene = rma.gene AND
        rm.position = rma.position AND
        rm.amino_acid = rma.amino_acid
    );

DELETE FROM _drm_articles a
  WHERE
    NOT EXISTS (
      SELECT 1 FROM resistance_mutations rm
      WHERE
        rm.gene = a.gene AND
        rm.position = a.position AND
        rm.amino_acid = a.amino_acid
    );

INSERT INTO resistance_mutation_articles
  SELECT gene, ref_name FROM _drm_articles
  GROUP BY gene, ref_name;

DROP TABLE approved_drugs;
DROP TABLE _drm_articles;
