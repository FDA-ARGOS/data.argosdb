CREATE VIEW IF NOT EXISTS mab_target_view
AS
SELECT
    DISTINCT *
FROM
    antibody_targets
WHERE
    ab_name NOT IN (
        SELECT
            ab_name
        FROM
            antibody_targets
        WHERE
            pdb_id IS NOT null
        )
UNION
SELECT
    DISTINCT *
FROM
    antibody_targets
WHERE
    pdb_id IS NOT null
;

CREATE VIEW IF NOT EXISTS mab_epitope_view
AS
SELECT
    ab_name,
    GROUP_CONCAT(position, '+') epitope
FROM
    antibody_epitopes
GROUP BY
    ab_name
;

CREATE VIEW IF NOT EXISTS mab_synonym
AS
SELECT
    ab_name,
    GROUP_CONCAT(synonym, ';') synonyms
FROM
    antibody_synonyms
GROUP BY
    ab_name
;

CREATE VIEW IF NOT EXISTS single_mab_view
AS
SELECT
    a.ab_name,
    s.synonyms,
    a.availability,
    t.pdb_id,
    t.target,
    t.class,
    e.epitope,
    a.institute,
    a.origin
FROM
    antibodies a
    LEFT JOIN
    mab_target_view t
    ON
    a.ab_name = t.ab_name
    LEFT JOIN
    mab_epitope_view e
    ON
    a.ab_name = e.ab_name
    LEFT JOIN
    mab_synonym s
    ON
    a.ab_name = s.ab_name
;
