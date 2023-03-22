-- for debuging functions, see https://dba.stackexchange.com/a/23357
-- LOAD 'auto_explain';
-- SET auto_explain.log_min_duration = 15000;
-- SET auto_explain.log_nested_statements = ON;


CREATE FUNCTION get_mutation_display(_mutobjs mutation_type[]) RETURNS VARCHAR[] AS $$
DECLARE
  _prev_del_gene VARCHAR;
  _prev_del_pos INTEGER;
  _prev_del_pos_start INTEGER;
  _prev_del VARCHAR;
  _mut VARCHAR;
  _gene VARCHAR;
  _pos INTEGER;
  _aa amino_acid_enum;
  _gene_display VARCHAR;
  _ref_aa amino_acid_enum;
  _mutations VARCHAR[];
BEGIN
  FOR _gene, _pos, _aa, _gene_display, _ref_aa IN
    SELECT mut.gene, mut.pos, mut.aa, g.display_name, ref.amino_acid
    FROM UNNEST(_mutobjs) mut
      JOIN genes AS g ON g.gene = mut.gene
      JOIN ref_amino_acid AS ref ON
        ref.gene = mut.gene AND
        ref.position = mut.pos
    ORDER BY g.gene_order, mut.pos
  LOOP
    IF _aa = 'del' THEN
      IF _prev_del_gene = _gene AND _prev_del_pos + 1 = _pos THEN
        _mut := 'Δ' || _prev_del_pos_start::varchar || '-' || _pos::varchar;
        _prev_del := NULL;
      ELSE
        _mut := 'Δ' || _pos::varchar;
        _prev_del_pos_start = _pos;
      END IF;
      _prev_del_gene = _gene;
      _prev_del_pos = _pos;
    ELSE
      _mut := _ref_aa::varchar || _pos::varchar || _aa::varchar;
    END IF;
    IF _gene != 'S' THEN
      _mut := _gene_display || ':' || _mut;
    END IF;
    IF _prev_del IS NOT NULL THEN
      _mutations := array_append(_mutations, _prev_del);
      _prev_del := NULL;
    END IF;
    IF _aa = 'del' THEN
      _prev_del := _mut;
    ELSE
      _mutations := array_append(_mutations, _mut);
    END IF;
  END LOOP;
  IF _prev_del IS NOT NULL THEN
    _mutations := array_append(_mutations, _prev_del);
  END IF;
  RETURN _mutations;
END
$$ LANGUAGE PLPGSQL IMMUTABLE;


CREATE FUNCTION get_indiv_mut_position(_iso_aggkey VARCHAR) RETURNS VARCHAR AS $$
DECLARE
  _mutobjs mutation_type[];
  _filtered mutation_type[];
BEGIN
  _mutobjs := get_isolate_agg_mutobjs(_iso_aggkey);
  _filtered := (
    SELECT ARRAY_AGG(mut)
    FROM UNNEST(_mutobjs) mut
    WHERE
      NOT EXISTS (
        SELECT 1 FROM ignore_mutations igm
        WHERE
          igm.gene = mut.gene AND
          igm.position = mut.pos AND
          igm.amino_acid = mut.aa
      )
  );
  IF ARRAY_LENGTH(_filtered, 1) = 0 THEN
    _filtered := _mutobjs;
  END IF;
  IF ARRAY_LENGTH(_filtered, 1) != 1 THEN
    RETURN NULL;
  ELSIF EXISTS (
    SELECT 1 FROM UNNEST(_filtered) WHERE aa = 'del'
  ) THEN
    RETURN NULL;
  ELSE
    RETURN _filtered[1].gene || ':' || _filtered[1].pos;
  END IF;
END
$$ LANGUAGE PLPGSQL IMMUTABLE;


CREATE FUNCTION get_isolate_display(_iso_name VARCHAR) RETURNS VARCHAR AS $$
DECLARE
  _mutobjs mutation_type[];
  _mutations VARCHAR[];
BEGIN
  _mutobjs := get_isolate_mutobjs('S', _iso_name);
  _mutations := get_mutation_display(_mutobjs);
  IF _mutations IS NULL OR array_length(_mutations, 1) = 0 THEN
    RETURN 'Wildtype';
  ELSE
    RETURN ARRAY_TO_STRING(_mutations, ' + ');
  END IF;
END
$$ LANGUAGE PLPGSQL IMMUTABLE;


SELECT
  iso_name,
  get_isolate_display(iso_name) AS iso_display
INTO TABLE isolate_displays
FROM isolates;

CREATE INDEX ON isolate_displays (iso_name);


CREATE FUNCTION get_isolate_agg_var_name(_iso_aggkey VARCHAR) RETURNS VARCHAR AS $$
  SELECT var_name
  FROM isolates iso
  WHERE var_name IS NOT NULL AND EXISTS(
    SELECT 1 FROM isolate_pairs pair
    WHERE pair.gene = 'S' AND pair.iso_name = iso.iso_name AND pair.iso_aggkey = _iso_aggkey
  )
  LIMIT 1
$$ LANGUAGE SQL IMMUTABLE;


CREATE FUNCTION get_isolate_agg_mutobjs(_iso_aggkey VARCHAR) RETURNS mutation_type[] AS $$
  SELECT ARRAY_AGG(
    DISTINCT (
      mut.gene,
      mut.position,
      mut.amino_acid
    )::mutation_type
  )
  FROM
    isolate_mutations mut,
    isolate_pairs pair
  WHERE
    mut.gene = 'S' AND
    pair.gene = 'S' AND
    pair.iso_aggkey = _iso_aggkey AND
    mut.iso_name = pair.iso_name AND
    /*NOT EXISTS (
      SELECT 1 FROM isolate_mutations ctl_mut
      WHERE
        ctl_mut.iso_name = pair.control_iso_name AND
        ctl_mut.gene = mut.gene AND
        ctl_mut.position = mut.position AND
        ctl_mut.amino_acid = mut.amino_acid
    ) AND*/
    NOT EXISTS (
      SELECT 1 FROM ignore_mutations igm
      WHERE
        igm.gene = mut.gene AND
        igm.position = mut.position AND
        igm.amino_acid = mut.amino_acid
    )
$$ LANGUAGE SQL IMMUTABLE;

CREATE FUNCTION get_isolate_agg_display(_iso_aggkey VARCHAR) RETURNS VARCHAR AS $$
DECLARE
  _mutobjs mutation_type[];
  _mutations VARCHAR[];
BEGIN
  _mutobjs := get_isolate_agg_mutobjs(_iso_aggkey);
  _mutations := get_mutation_display(_mutobjs);
  IF _mutations IS NULL OR ARRAY_LENGTH(_mutations, 1) = 0 THEN
    RETURN 'Wildtype';
  ELSE
    RETURN ARRAY_TO_STRING(_mutations, ' + ');
  END IF;
END
$$ LANGUAGE PLPGSQL IMMUTABLE;

SELECT DISTINCT
  iso_aggkey,
  ARRAY_LENGTH(
    get_countable_mutations(
      get_isolate_agg_mutobjs(iso_aggkey)
    ),
    1
  ) AS num_mutations,
  get_indiv_mut_position(iso_aggkey) AS position,
  get_isolate_agg_display(iso_aggkey) AS iso_agg_display,
  get_isolate_agg_var_name(iso_aggkey) AS var_name
INTO TABLE isolate_aggkeys
FROM isolate_pairs
WHERE gene = 'S';

CREATE INDEX ON isolate_aggkeys (iso_aggkey);

SELECT
  ref_name,
  rx_name,
  csv_agg(ab.ab_name ORDER BY ab.priority)::VARCHAR AS antibody_names,
  BOOL_AND(ab.visibility) AS visibility,
  MAX(ab.priority) * 10 + COUNT(ab.ab_name) AS antibody_order
INTO TABLE rx_antibody_names_all
FROM rx_antibodies rxab, antibodies ab
WHERE
  rxab.ab_name = ab.ab_name
GROUP BY ref_name, rx_name;

CREATE INDEX ON rx_antibody_names_all (ref_name, rx_name);

SELECT
  ref_name,
  rx_name,
  csv_agg(ab.ab_name ORDER BY ab.priority)::VARCHAR AS antibody_names,
  BOOL_AND(ab.visibility) AS visibility,
  MAX(ab.priority) * 10 + COUNT(ab.ab_name) AS antibody_order,
  COUNT(ab.ab_name) AS num_antibodies
INTO TABLE rx_antibody_names_any
FROM rx_antibodies rxab, antibodies ab
WHERE
  rxab.ab_name = ab.ab_name
GROUP BY ref_name, rx_name;

ALTER TABLE rx_antibody_names_any
ADD CONSTRAINT PK_rx_antibody_names_any PRIMARY KEY (ref_name, rx_name, antibody_names);

SELECT
  DISTINCT
  antibody_names,
  ab_name,
  visibility,
  antibody_order,
  num_antibodies
INTO TABLE tmp_antibody_names
FROM rx_antibody_names_any, rx_antibodies rxab
WHERE
  rx_antibody_names_any.ref_name = rxab.ref_name AND
  rx_antibody_names_any.rx_name = rxab.rx_name AND
  num_antibodies > 1;

INSERT INTO rx_antibody_names_any (
  ref_name,
  rx_name,
  antibody_names,
  visibility,
  antibody_order,
  num_antibodies
) SELECT
  DISTINCT
  ref_name,
  rx_name,
  antibody_names,
  visibility,
  antibody_order,
  num_antibodies
FROM rx_antibodies rxab, tmp_antibody_names abns
WHERE
  rxab.ab_name = abns.ab_name
ON CONFLICT DO NOTHING;


CREATE TYPE susc_summary_agg_key AS ENUM (
  'rx_type',
  'article',
  'infected_variant',
  'vaccine',
  'antibody',
  'antibody:any',
  'antibody:indiv',
  'variant',
  'isolate_agg',
  'isolate',
  'position',
  'vaccine_dosage',
  'timing',
  'subject_species',
  'potency_type',
  'potency_unit'
);


CREATE FUNCTION summarize_susc_results(_agg_by susc_summary_agg_key[]) RETURNS VOID AS $$
  DECLARE
    _stmt TEXT;
    _plasma_related_agg_by susc_summary_agg_key[];
    _isolate_related_agg_by susc_summary_agg_key[];
    _ext_col_names TEXT[];
    _ext_col_values TEXT[];
    _ext_joins TEXT[];
    _ext_where TEXT[];
    _ext_group_by TEXT[];
    _ret INTEGER;
  BEGIN
    _stmt := $STMT$
      WITH rows AS(
        INSERT INTO susc_summary (
          aggregate_by,
          num_studies,
          num_samples,
          num_experiments
          %1s
        ) SELECT
          $1 AS aggregate_by,
          COUNT(DISTINCT S.ref_name) AS num_studies,
          unique_sum(
            ARRAY_AGG((
              S.ref_name || '$##$' || S.rx_group,
              S.cumulative_count
            )::unique_sum_type)
          ) AS num_samples,
          unique_sum(
            ARRAY_AGG((
              -- XXX: This uniqkey list must be identical to the primary keys
              -- of table susc_results
              S.ref_name || '$##$' ||
              S.rx_group || '$##$' ||
              S.control_iso_name || '$##$' ||
              S.iso_name || '$##$' ||
              S.potency_type || '$##$' ||
              S.control_assay_name || '$##$' ||
              S.assay_name,
              S.cumulative_count
            )::unique_sum_type)
          ) AS num_experiments
          %2s
        FROM
          susc_results S %3s
        WHERE %4s
        GROUP BY aggregate_by %5s
        RETURNING 1
      )
      SELECT count(*) FROM rows;
    $STMT$;

    IF 'rx_type' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'rx_type');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'S.rx_type');
      _ext_group_by := ARRAY_APPEND(_ext_group_by, 'S.rx_type');
    END IF;

    IF 'article' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'ref_name');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'S.ref_name');
      _ext_group_by := ARRAY_APPEND(_ext_group_by, 'S.ref_name');
    END IF;

    IF 'antibody:indiv' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        rx_type,
        antibody_names,
        antibody_order
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        'antibody' AS rx_type,
        rxab.ab_name AS antibody_names,
        ab.priority AS antibody_order
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN rx_antibodies rxab ON
          S.ref_name = rxab.ref_name AND
          S.rx_name = rxab.rx_name
        JOIN antibodies ab ON rxab.ab_name = ab.ab_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        rxab.ab_name,
        ab.priority
      $X$);
    END IF;

    IF 'antibody:any' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        rx_type,
        antibody_names,
        antibody_order
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        'antibody' AS rx_type,
        ab.antibody_names AS antibody_names,
        ab.antibody_order AS antibody_order
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN rx_antibody_names_any ab ON
          S.ref_name = ab.ref_name AND
          S.rx_name = ab.rx_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        ab.antibody_names,
        ab.antibody_order
      $X$);
    END IF;

    IF 'antibody' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        rx_type,
        antibody_names,
        antibody_order
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        'antibody' AS rx_type,
        ab.antibody_names AS antibody_names,
        ab.antibody_order AS antibody_order
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN rx_antibody_names_all ab ON
          S.ref_name = ab.ref_name AND
          S.rx_name = ab.rx_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        ab.antibody_names,
        ab.antibody_order
      $X$);
    END IF;

    IF 'vaccine' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        rx_type,
        vaccine_name,
        vaccine_order,
        num_subjects
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        'vacc-plasma' AS rx_type,
        rxvp.vaccine_name,
        v.priority AS vaccine_order,
        unique_sum(
          ARRAY_AGG((
            sbj.ref_name || '$##$' || sbj.subject_name,
            sbj.num_subjects
          )::unique_sum_type)
        ) AS num_subjects
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN rx_vacc_plasma rxvp ON
          S.ref_name = rxvp.ref_name AND (
            S.rx_name = rxvp.rx_name OR
            EXISTS (
              SELECT 1 FROM unlinked_susc_results usr
              WHERE
                S.ref_name = usr.ref_name AND
                S.rx_group = usr.rx_group AND
                usr.rx_name = rxvp.rx_name
            )
          )
        JOIN vaccines v ON
          rxvp.vaccine_name = v.vaccine_name
        JOIN subjects sbj ON
          S.ref_name = sbj.ref_name AND
          rxvp.subject_name = sbj.subject_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        rxvp.vaccine_name,
        v.priority
      $X$);
    END IF;

    IF 'vaccine_dosage' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'vaccine_dosage');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'rxvp.dosage');
      _ext_group_by := ARRAY_APPEND(_ext_group_by, 'rxvp.dosage');
    END IF;

    IF 'infected_variant' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        rx_type,
        infected_var_name,
        num_subjects
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        'conv-plasma' AS rx_type,
        CASE
          WHEN infected.as_wildtype IS TRUE THEN 'Wild Type'
          ELSE infected.var_name
        END AS infected_var_name,
        unique_sum(
          ARRAY_AGG((
            sbj.ref_name || '$##$' || sbj.subject_name,
            sbj.num_subjects
          )::unique_sum_type)
        ) AS num_subjects
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN rx_conv_plasma rx ON
          S.ref_name = rx.ref_name AND (
            S.rx_name = rx.rx_name OR
            EXISTS (
              SELECT 1 FROM unlinked_susc_results usr
              WHERE
                S.ref_name = usr.ref_name AND
                S.rx_group = usr.rx_group AND
                usr.rx_name = rx.rx_name
            )
          )
        JOIN subjects sbj ON
          rx.ref_name = sbj.ref_name AND
          rx.subject_name = sbj.subject_name
        LEFT JOIN variants infected ON
          rx.infected_var_name = infected.var_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        infected.as_wildtype,
        infected.var_name
      $X$);
    END IF;

    IF 'timing' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'timing');
      IF 'vaccine' = ANY(_agg_by) THEN
        _ext_col_values := ARRAY_APPEND(_ext_col_values, 'rxvp.timing');
        _ext_group_by := ARRAY_APPEND(_ext_group_by, 'rxvp.timing');
      ELSE
        _ext_col_values := ARRAY_APPEND(_ext_col_values, 'rx.timing');
        _ext_group_by := ARRAY_APPEND(_ext_group_by, 'rx.timing');
      END IF;
    END IF;

    _plasma_related_agg_by := ARRAY['vaccine', 'infected_variant'];
    _isolate_related_agg_by := ARRAY['variant', 'isolate_agg', 'isolate'];

    IF _isolate_related_agg_by && _agg_by THEN
      -- query num_subjects when any _plasma_related_agg_by is not used
      IF NOT (_plasma_related_agg_by && _agg_by) THEN
        _ext_col_names := ARRAY_APPEND(_ext_col_names, 'num_subjects');
        _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
          unique_sum(
            ARRAY_AGG((
              sbj.ref_name || '$##$' || sbj.subject_name,
              sbj.num_subjects
            )::unique_sum_type)
          ) AS num_subjects
        $X$);
        _ext_joins := ARRAY_APPEND(_ext_joins, $X$
          LEFT JOIN subject_plasma sbjp ON
            S.ref_name = sbjp.ref_name AND (
              S.rx_name = sbjp.rx_name OR
              EXISTS (
                SELECT 1 FROM unlinked_susc_results usr
                WHERE
                  S.ref_name = usr.ref_name AND
                  S.rx_group = usr.rx_group AND
                  usr.rx_name = sbjp.rx_name
              )
            )
          LEFT JOIN subjects sbj ON
            sbjp.ref_name = sbj.ref_name AND
            sbjp.subject_name = sbj.subject_name
        $X$);
      END IF;
    END IF;

    IF 'variant' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'var_name');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'target.var_name');
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN isolates target ON
          S.iso_name = target.iso_name
      $X$);
      _ext_where := ARRAY_APPEND(_ext_where, $X$
        target.var_name IS NOT NULL
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, 'target.var_name');
    END IF;

    IF 'isolate_agg' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        iso_type,
        iso_aggkey,
        iso_agg_display,
        var_name
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        CASE WHEN isoagg.num_mutations = 1 THEN
          'indiv-mut'
        ELSE
          'combo-muts'
        END::iso_type_enum AS iso_type,
        isoagg.iso_aggkey AS iso_aggkey,
        isoagg.iso_agg_display AS iso_agg_display,
        isoagg.var_name AS var_name
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN isolate_pairs pair ON
          pair.gene = 'S' AND
          S.control_iso_name = pair.control_iso_name AND
          S.iso_name = pair.iso_name
        JOIN isolate_aggkeys isoagg ON
          pair.iso_aggkey = isoagg.iso_aggkey
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        isoagg.num_mutations,
        isoagg.iso_aggkey,
        isoagg.iso_agg_display,
        isoagg.var_name
      $X$);
    END IF;

    IF 'potency_type' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        potency_type
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        S.potency_type AS potency_type
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        S.potency_type
      $X$);
    END IF;

    IF 'potency_unit' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        potency_unit
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        S.potency_unit AS potency_unit
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        S.potency_unit
      $X$);
    END IF;

    IF 'isolate' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, $X$
        iso_type,
        control_iso_name,
        control_iso_display,
        control_var_name,
        iso_name,
        iso_display,
        var_name,
        iso_aggkey,
        iso_agg_display
      $X$);
      _ext_col_values := ARRAY_APPEND(_ext_col_values, $X$
        CASE WHEN pair.num_mutations = 1 THEN
          'indiv-mut'
        ELSE
          'combo-muts'
        END::iso_type_enum AS iso_type,
        pair.control_iso_name,
        ctl_display.iso_display,
        control.var_name,
        pair.iso_name,
        tgt_display.iso_display,
        target.var_name,
        pair.iso_aggkey,
        isoagg.iso_agg_display
      $X$);
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN isolate_pairs pair ON
          pair.gene = 'S' AND
          S.iso_name = pair.iso_name AND
          S.control_iso_name = pair.control_iso_name
        LEFT JOIN isolate_aggkeys isoagg ON
          pair.iso_aggkey = isoagg.iso_aggkey
        JOIN isolates target ON
          S.iso_name = target.iso_name
        JOIN isolates control ON
          S.control_iso_name = control.iso_name
        JOIN isolate_displays tgt_display ON
          S.iso_name = tgt_display.iso_name
        JOIN isolate_displays ctl_display ON
          S.control_iso_name = ctl_display.iso_name
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        pair.num_mutations,
        pair.control_iso_name,
        ctl_display.iso_display,
        control.var_name,
        pair.iso_name,
        tgt_display.iso_display,
        target.var_name,
        pair.iso_aggkey,
        isoagg.iso_agg_display
      $X$);
    END IF;

    IF 'position' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'position');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'isoagg.position');
      _ext_joins := ARRAY_APPEND(_ext_joins, $X$
        JOIN isolate_pairs pair ON
          pair.gene = 'S' AND
          S.control_iso_name = pair.control_iso_name AND
          S.iso_name = pair.iso_name
        JOIN isolate_aggkeys isoagg ON
          pair.iso_aggkey = isoagg.iso_aggkey
      $X$);
      _ext_where := ARRAY_APPEND(_ext_where, $X$
        isoagg.position IS NOT NULL
      $X$);
      _ext_group_by := ARRAY_APPEND(_ext_group_by, $X$
        isoagg.position
      $X$);
    END IF;

    IF 'subject_species' = ANY(_agg_by) THEN
      _ext_col_names := ARRAY_APPEND(_ext_col_names, 'subject_species');
      _ext_col_values := ARRAY_APPEND(_ext_col_values, 'subject_species');
      _ext_group_by := ARRAY_APPEND(_ext_group_by, 'subject_species');
    END IF;

    IF _ext_where IS NULL THEN
      _ext_where := ARRAY['TRUE'];
    END IF;
    
    _stmt := FORMAT(
      _stmt,
      CASE WHEN _ext_col_names IS NULL THEN ''
      ELSE ', ' || ARRAY_TO_STRING(_ext_col_names, ', ') END,
      CASE WHEN _ext_col_values IS NULL THEN ''
      ELSE ', ' || ARRAY_TO_STRING(_ext_col_values, ', ') END,
      ARRAY_TO_STRING(_ext_joins, ' '),
      ARRAY_TO_STRING(_ext_where, ' AND '),
      CASE WHEN _ext_group_by IS NULL THEN ''
      ELSE ', ' || ARRAY_TO_STRING(_ext_group_by, ', ') END
    );

    EXECUTE _stmt INTO _ret USING ARRAY_TO_STRING((
      SELECT ARRAY_AGG(x) FROM (
        SELECT UNNEST(_agg_by) AS x
        ORDER BY x
      ) AS _
    ), ',');
    RAISE NOTICE
      'Summarized susc_results by %: % rows created',
      ARRAY_TO_STRING(_agg_by, ' + '),
      _ret;
  END
$$ LANGUAGE PLPGSQL VOLATILE;


DO $$
  DECLARE
    _agg_by_auto_options susc_summary_agg_key[];
    _agg_by susc_summary_agg_key[];
    _rx_agg_by susc_summary_agg_key[];
    _iso_agg_by susc_summary_agg_key[];
  BEGIN
    _agg_by_auto_options := ARRAY[
      'rx_type',
      'article',
      'infected_variant',
      'vaccine',
      'antibody',
      'antibody:any',
      'antibody:indiv',
      'variant',
      'isolate_agg',
      'isolate',
      'position'
    ];
    _rx_agg_by := ARRAY[
      'rx_type',
      'antibody',
      'antibody:any',
      'antibody:indiv',
      'vaccine',
      'infected_variant'
    ];
    _iso_agg_by := ARRAY[
      'isolate_agg',
      'variant',
      'isolate',
      'position'
    ];

    FOR _agg_by IN
    SELECT agg_by.agg_by FROM (

      -- combination of one element
      SELECT DISTINCT ARRAY[one] agg_by
      FROM UNNEST(_agg_by_auto_options) one

      UNION
      -- combination of two elements
      SELECT DISTINCT ARRAY[one, two] agg_by
      FROM
        UNNEST(_agg_by_auto_options) one,
        UNNEST(_agg_by_auto_options) two
      WHERE
        one < two AND
        array_intersect_count(ARRAY[one, two], _rx_agg_by) < 2 AND
        array_intersect_count(ARRAY[one, two], _iso_agg_by) < 2

      UNION
      -- combination of three elements
      SELECT DISTINCT ARRAY[one, two, three] agg_by
      FROM
        UNNEST(_agg_by_auto_options) one,
        UNNEST(_agg_by_auto_options) two,
        UNNEST(_agg_by_auto_options) three
      WHERE
        one < two AND two < three AND
        array_intersect_count(ARRAY[one, two, three], _rx_agg_by) < 2 AND
        array_intersect_count(ARRAY[one, two, three], _iso_agg_by) < 2
    ) agg_by
    ORDER BY agg_by.agg_by
    LOOP
      PERFORM summarize_susc_results(_agg_by);
    END LOOP;

    /*PERFORM summarize_susc_results(ARRAY[
      'article',
      'infected_variant',
      'isolate',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'article',
      'infected_variant',
      'isolate',
      'timing',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'article',
      'vaccine',
      'isolate',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'article',
      'vaccine',
      'isolate',
      'timing',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'antibody',
      'variant',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'antibody',
      'isolate_agg',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'infected_variant',
      'variant',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'infected_variant',
      'variant',
      'timing',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'infected_variant',
      'isolate_agg',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'infected_variant',
      'isolate_agg',
      'timing',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'vaccine',
      'variant',
      'vaccine_dosage',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'vaccine',
      'variant',
      'vaccine_dosage',
      'timing',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'vaccine',
      'isolate_agg',
      'vaccine_dosage',
      'subject_species',
      'potency_type',
      'potency_unit'
    ]::susc_summary_agg_key[]);

    PERFORM summarize_susc_results(ARRAY[
      'vaccine',
      'isolate_agg',
      'vaccine_dosage',
      'timing',
      'subject_species',
      'potency_unit',
      'potency_type'
    ]::susc_summary_agg_key[]);*/

    PERFORM summarize_susc_results(ARRAY[]::susc_summary_agg_key[]);
  END
$$ LANGUAGE PLPGSQL;

DROP TABLE isolate_displays;
DROP TABLE isolate_aggkeys;
DROP TABLE tmp_antibody_names;
DROP TABLE rx_antibody_names_any;
DROP TABLE rx_antibody_names_all;
DROP FUNCTION get_isolate_mutobjs;
DROP FUNCTION get_indiv_mut_position;
DROP FUNCTION get_isolate_agg_var_name;
DROP FUNCTION get_isolate_agg_mutobjs;
DROP FUNCTION get_isolate_agg_display;
DROP FUNCTION summarize_susc_results;
