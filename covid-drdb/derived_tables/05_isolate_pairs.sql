-- for debuging functions, see https://dba.stackexchange.com/a/23357
-- LOAD 'auto_explain';
-- SET auto_explain.log_min_duration = 15000;
-- SET auto_explain.log_nested_statements = ON;


CREATE TYPE mutation_type AS (
  gene VARCHAR,
  pos INTEGER,
  aa amino_acid_enum
);


CREATE FUNCTION get_countable_mutations(_mutobjs mutation_type[]) RETURNS VARCHAR[] AS $$
DECLARE
  _prev_gene varchar;
  _gene varchar;
  _gene_order integer;
  _pos integer;
  _pos_end integer;
  _aa amino_acid_enum;
  _mut varchar;
  _mutations varchar[];
BEGIN
  FOR _gene, _gene_order, _pos, _pos_end, _aa IN
    SELECT
      DISTINCT
      mut.gene,
      g.gene_order,
      CASE WHEN
        mut.aa = 'del' AND
        position_start IS NOT NULL
      THEN
        position_start
      ELSE
        mut.pos
      END AS pos_start,
      position_end AS pos_end,
      mut.aa
    FROM UNNEST(_mutobjs) mut
      JOIN genes AS g ON g.gene = mut.gene
      LEFT JOIN known_deletion_ranges dr ON
        dr.gene = mut.gene AND
        mut.pos BETWEEN position_start AND position_end
    WHERE
      NOT EXISTS (
        SELECT 1 FROM ignore_mutations igm
        WHERE
          igm.gene = mut.gene AND
          igm.position = mut.pos AND
          igm.amino_acid = mut.aa
      )
    ORDER BY g.gene_order, pos_start
  LOOP
    IF _aa = 'del' AND _pos_end IS NOT NULL THEN
      _mut := _pos::VARCHAR || '-' || _pos_end::VARCHAR || _aa::VARCHAR;
    ELSE
      _mut := _pos::VARCHAR || _aa::VARCHAR;
    END IF;
    IF _prev_gene IS NULL OR _prev_gene != _gene THEN
      _mut := _gene || ':' || _mut;
    END IF;
    _prev_gene := _gene;
    _mutations := array_append(_mutations, _mut);
  END LOOP;
  RETURN _mutations;
END
$$ LANGUAGE PLPGSQL IMMUTABLE;


CREATE FUNCTION get_isolate_aggkey(_gene VARCHAR, my_iso_name VARCHAR, my_control_iso_name VARCHAR) RETURNS VARCHAR AS $$
DECLARE
  _mutobjs mutation_type[];
  _mutations varchar[];
BEGIN
  IF my_control_iso_name IS NULL THEN
    my_control_iso_name := 'Wuhan-Hu-1';
  END IF;
  _mutobjs := (
    SELECT ARRAY_AGG(
      (
        mut.gene,
        mut.position,
        mut.amino_acid
      )::mutation_type
    )
    FROM isolate_mutations mut
    WHERE
      mut.gene = _gene AND
      mut.iso_name = my_iso_name/* AND
      NOT EXISTS (
        SELECT 1 FROM isolate_mutations ctl_mut
        WHERE
          ctl_mut.iso_name = my_control_iso_name AND
          ctl_mut.gene = mut.gene AND
          ctl_mut.position = mut.position AND
          ctl_mut.amino_acid = mut.amino_acid
      ) */
  );
  _mutations := get_countable_mutations(_mutobjs);
  RETURN ARRAY_TO_STRING(_mutations, '+');
END
$$ LANGUAGE PLPGSQL IMMUTABLE;

CREATE FUNCTION get_isolate_mutobjs(_gene VARCHAR, _iso_name VARCHAR) RETURNS mutation_type[] AS $$
  SELECT ARRAY_AGG(
    (
      mut.gene,
      mut.position,
      mut.amino_acid
    )::mutation_type
  )
  FROM isolate_mutations mut
  WHERE
    mut.gene = _gene AND
    mut.iso_name = _iso_name
$$ LANGUAGE SQL IMMUTABLE;

CREATE FUNCTION populate_isolate_pairs(_gene VARCHAR) RETURNS VOID AS $$
INSERT INTO isolate_pairs (
  gene,
  control_iso_name,
  iso_name,
  iso_aggkey,
  num_mutations
) SELECT DISTINCT
  _gene AS gene,
  control_iso_name,
  iso_name,
  get_isolate_aggkey(_gene, iso_name, control_iso_name) AS iso_aggkey,
  ARRAY_LENGTH(
    get_countable_mutations(
      get_isolate_mutobjs(_gene, iso_name)
    ),
    1
  ) AS num_mutations
FROM susc_results;
$$ LANGUAGE SQL VOLATILE;

DO $$
BEGIN
  PERFORM populate_isolate_pairs('_3CLpro');
  PERFORM populate_isolate_pairs('RdRP');
  PERFORM populate_isolate_pairs('S');
END
$$ LANGUAGE PLPGSQL;

DROP FUNCTION get_isolate_aggkey;
DROP FUNCTION populate_isolate_pairs;
