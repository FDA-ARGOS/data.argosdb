CREATE VIEW IF NOT EXISTS isolate_mutations_with_ref_view
AS
SELECT
    a.iso_name,
    a.gene,
    b.amino_acid as ref,
    a.position,
    CASE
        WHEN a.amino_acid == 'del' THEN
            '∆'
        ELSE
            a.amino_acid
    END amino_acid,
    CASE
        WHEN a.gene == 'S' AND a.position >= 1 AND a.position <= 305 THEN
            'NTD'
        WHEN a.gene == 'S' AND a.position >= 306 AND a.position <= 534 THEN
            'RBD'
        WHEN a.gene == 'S' AND a.position >= 535 AND a.position <= 686 THEN
            'CTD'
        WHEN a.gene == 'S' AND a.position >= 687 AND a.position <= 1273 THEN
            'S2'
        ELSE
            ''
    END domain,
    b.amino_acid || a.position || a.amino_acid single_mut_name
FROM
    isolate_mutations a,
    ref_amino_acid b
WHERE
    a.gene = b.gene
    AND
    a.position = b.position
;

CREATE VIEW IF NOT EXISTS isolate_mutations_spike_view
AS
SELECT
    *
FROM
    isolate_mutations_with_ref_view
WHERE
    gene == 'S'
;

CREATE VIEW IF NOT EXISTS isolate_mutations_spike_no_614g_683g_view
AS
SELECT
    *
FROM
    isolate_mutations_spike_view
WHERE
    (position, amino_acid) != (614, 'G')
    AND
    (position, amino_acid) != (683, 'G')
;

-- merge NTD mutations
CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_del_view
AS
SELECT
    *
FROM
    isolate_mutations_spike_no_614g_683g_view
WHERE
    amino_acid = '∆'
;

CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_del_group_view
AS
WITH firstrows AS
(
    SELECT
        a.iso_name,
        a.position,
        ROW_NUMBER() OVER (PARTITION BY iso_name ORDER BY position) group_id
    FROM
        _isolate_mutations_spike_del_view a
    WHERE
        NOT EXISTS (
            SELECT 1
            FROM _isolate_mutations_spike_del_view b
            WHERE
                a.iso_name = b.iso_name
                AND
                (b.position - a.position) = 1
        )

)
SELECT *,
    (
        SELECT b.group_id
        FROM firstrows b
        WHERE
            a.position <= b.position
            AND
            b.iso_name = a.iso_name
    ) group_id
FROM
_isolate_mutations_spike_del_view a
;


CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_consecu_del_view
AS
SELECT
    a.iso_name,
    min(a.position) min_pos,
    max(a.position) max_pos,
    CASE
        WHEN min(a.position) == max(a.position) THEN
            min(a.position)
        ELSE
            min(a.position) || '-' || max(a.position)
    END pos_range,
    (max(a.position) - min(a.position)) pos_range_length
FROM
    _isolate_mutations_spike_del_group_view a
GROUP BY
    a.iso_name,
    a.group_id
ORDER BY
    a.iso_name,
    a.position
;


CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_merged_del_view
AS
SELECT
    a.iso_name,
    a.gene,
    group_concat(a.ref, '') ref,
    b.min_pos position,
    -- b.min_pos start_pos,
    -- b.max_pos stop_pos,
    a.amino_acid,
    a.domain domain,
    group_concat(a.ref, '') || b.pos_range || '∆' single_mut_name
from
    _isolate_mutations_spike_del_view a,
    _isolate_mutations_spike_consecu_del_view b
WHERE
    a.iso_name = b.iso_name
    AND
    a.position >= b.min_pos
    AND
    a.position <= b.max_pos
group by
    a.iso_name,
    b.pos_range

UNION

SELECT
    a.*
from
    _isolate_mutations_spike_del_view a
WHERE
    NOT EXISTS
    (
        SELECT 1
        FROM _isolate_mutations_spike_consecu_del_view b
        WHERE
            a.iso_name = b.iso_name
            AND
            a.position >= b.min_pos
            AND
            a.position <= b.max_pos
    )
;

CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_normal_view
AS
SELECT
    *
FROM
    _isolate_mutations_spike_merged_del_view

UNION

SELECT
    *
FROM
    isolate_mutations_spike_no_614g_683g_view
WHERE
    amino_acid != '∆'
;


CREATE VIEW IF NOT EXISTS isolate_mutations_single_s_mut_view
AS
SELECT
    *
FROM
    _isolate_mutations_spike_normal_view a
WHERE
    EXISTS (
        SELECT
            1
        FROM
            _isolate_mutations_spike_normal_view b
        WHERE
            a.iso_name = b.iso_name
        GROUP BY
            b.iso_name
        HAVING
            COUNT(1) = 1
    )
;


CREATE VIEW IF NOT EXISTS _isolate_mutations_combo_s_mut_view
AS
SELECT
    *
FROM
    _isolate_mutations_spike_normal_view a
WHERE
    EXISTS (
        SELECT
            1
        FROM
            _isolate_mutations_spike_normal_view b
        WHERE
            a.iso_name = b.iso_name
        GROUP BY
            b.iso_name
        HAVING
            COUNT(1) > 1
    )
;


CREATE VIEW IF NOT EXISTS isolate_mutations_combo_s_mut_view
AS
SELECT
    iso.iso_name,
    iso.var_name,
    '' domain,
    GROUP_CONCAT(muts.single_mut_name, '+') pattern,
    COUNT(muts.single_mut_name) num_muts
FROM
    _isolate_mutations_combo_s_mut_view muts,
    isolates iso
WHERE
    muts.iso_name = iso.iso_name
GROUP BY
    iso.iso_name
ORDER BY
    iso.iso_name,
    muts.position
;


CREATE VIEW IF NOT EXISTS isolate_wildtype_view
AS
SELECT
    iso.*
FROM
    isolates iso,
    variants var
WHERE
    iso.var_name = var.var_name
    AND
    var.as_wildtype
;

CREATE VIEW IF NOT EXISTS isolate_non_wildtype_view
AS
SELECT
    iso.*
FROM
    isolates iso,
    variants var
WHERE
    iso.var_name = var.var_name
    AND
    NOT var.as_wildtype
UNION
SELECT
    iso.*
FROM
    isolates iso
WHERE
    iso.var_name IS NULL
;

CREATE VIEW IF NOT EXISTS isolate_omicron_view
AS
SELECT
    iso.*
FROM
    isolates iso,
    variants var
WHERE
    iso.var_name = var.var_name
    AND
    var.var_name LIKE 'Omicron/%'
;


CREATE VIEW IF NOT EXISTS isolate_omicron_ba2_view
AS
SELECT
    iso.*
FROM
    isolates iso,
    variants var
WHERE
    iso.var_name = var.var_name
    AND
    var.var_name = 'Omicron/BA.2'
;


CREATE VIEW IF NOT EXISTS isolate_variant_view
AS
SELECT
    iso.*,
    var.as_wildtype
FROM
    isolates iso
LEFT JOIN
    variants var
ON
    iso.var_name = var.var_name
;


CREATE VIEW IF NOT EXISTS _isolate_mutations_spike_raw_mut_view
AS
SELECT
    *
FROM
    _isolate_mutations_spike_merged_del_view

UNION

SELECT
    *
FROM
    isolate_mutations_spike_view
WHERE
    amino_acid != '∆'
;


CREATE VIEW IF NOT EXISTS _isolate_mutations_variant_raw_s_mut_view
AS
SELECT
    *
FROM
    _isolate_mutations_spike_raw_mut_view a
WHERE
    EXISTS (
        SELECT
            1
        FROM
            isolates b
        WHERE
            a.iso_name = b.iso_name
            AND
            b.var_name IS NOT NULL

    )
;


CREATE VIEW IF NOT EXISTS isolate_mutations_variant_raw_s_mut_view
AS
SELECT
    iso.iso_name,
    iso.var_name,
    '' domain,
    GROUP_CONCAT(muts.single_mut_name, '+') pattern
FROM
    _isolate_mutations_variant_raw_s_mut_view muts,
    isolates iso
WHERE
    muts.iso_name = iso.iso_name
GROUP BY
    iso.iso_name
ORDER BY
    iso.iso_name,
    muts.position
;
