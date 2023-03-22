SELECT
  *,
  -- when comparing two collection_date, we need to take into account the
  -- effect of collection_date_cmp. e.g. <2021-05-01 should smaller than
  -- =2021-05-01.
  CASE collection_date_cmp
    WHEN '<' THEN collection_date - 1
    WHEN '>' THEN collection_date + 1
    ELSE collection_date
  END comparable_collection_date
INTO dup_subject_isolates
FROM subject_isolates;

CREATE INDEX ON dup_subject_isolates (iso_name);
CREATE INDEX ON dup_subject_isolates (ref_name, subject_name);
CREATE INDEX ON dup_subject_isolates (comparable_collection_date);

INSERT INTO invivo_selection_results
  SELECT
    DISTINCT
    SbjIso.ref_name,
    SbjIso.subject_name,
    SbjInf.infected_var_name,
    IsoMuts.gene,
    IsoMuts.position,
    IsoMuts.amino_acid,
    CASE WHEN IsoMuts.count IS NULL THEN 1 ELSE IsoMuts.count END,
    CASE WHEN IsoMuts.total IS NULL THEN 1 ELSE IsoMuts.total END,
    SbjInf.infection_date,
    SbjIso.collection_date AS appearance_date,
    SbjInf.immune_status,
    SbjInf.severity
  FROM dup_subject_isolates SbjIso
  JOIN isolate_mutations IsoMuts ON
    SbjIso.iso_name = IsoMuts.iso_name
  JOIN subject_infections SbjInf ON
    SbjIso.ref_name = SbjInf.ref_name AND
    SbjIso.subject_name = SbjInf.subject_name
  WHERE
    EXISTS (
      SELECT 1
        FROM dup_subject_isolates SbjIsoPrev
        WHERE
          SbjIsoPrev.ref_name = SbjIso.ref_name AND
          SbjIsoPrev.subject_name = SbjIso.subject_name AND
          -- use comparable_collection_date to take into account of the collection_date_cmp
          SbjIsoPrev.comparable_collection_date < SbjIso.comparable_collection_date
    ) AND
    NOT EXISTS (
      SELECT 1
        FROM dup_subject_isolates SbjIsoPrev
        JOIN isolate_mutations IsoMutsPrev ON
          SbjIsoPrev.iso_name = IsoMutsPrev.iso_name
        WHERE
          SbjIsoPrev.ref_name = SbjIso.ref_name AND
          SbjIsoPrev.subject_name = SbjIso.subject_name AND
          -- use comparable_collection_date to take into account of the collection_date_cmp
          SbjIsoPrev.comparable_collection_date < SbjIso.comparable_collection_date AND
          IsoMutsPrev.gene = IsoMuts.gene AND
          IsoMutsPrev.position = IsoMuts.position AND
          IsoMutsPrev.amino_acid = IsoMuts.amino_acid
    );

DROP TABLE dup_subject_isolates;
