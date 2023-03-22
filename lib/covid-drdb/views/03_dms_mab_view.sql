CREATE VIEW IF NOT EXISTS dms_single_mab_view
AS
SELECT
    rx.ref_name ref_name,
    rx.rx_name rx_name,
    rx.ab_name ab_name,
    mab.synonyms synonyms,
    mab.availability,
    mab.pdb_id,
    mab.target,
    mab.class,
    mab.epitope,
    mab.institute,
    mab.origin
FROM
    rx_antibodies rx,
    single_mab_view mab
WHERE
    rx.ab_name = mab.ab_name
GROUP BY
    rx.ref_name,
    rx.rx_name
HAVING
    count(rx.ab_name) = 1
;

CREATE VIEW IF NOT EXISTS dms_combo_mab_view
AS
SELECT
    a.ref_name,
    a.rx_name,
    a.ab_name || '/' || b.ab_name ab_name,
    '' synonyms,
    CASE
        WHEN a.availability == b.availability THEN
            a.availability
        ELSE
            ''
    END availability,
    '' AS pdb_id,
    '' AS target,
    '' AS class,
    '' AS epitope,
    CASE
        WHEN a.institute == b.institute THEN
            a.institute
        ELSE
            ''
    END institute,
    CASE
        WHEN a.origin == b.origin THEN
            a.origin
        ELSE
            ''
    END origin
FROM
    (SELECT
        rx.*,
        ab.availability,
        ab.institute,
        ab.origin
    FROM
        rx_antibodies rx,
        antibodies ab
    WHERE
        rx.ab_name = ab.ab_name
    ) a,
    (SELECT
        rx.*,
        ab.availability,
        ab.institute,
        ab.origin
    FROM
        rx_antibodies rx,
        antibodies ab
    WHERE
        rx.ab_name = ab.ab_name
    ) b
WHERE
    a.ref_name = b.ref_name and
    a.rx_name = b.rx_name and
    a.ab_name != b.ab_name and
    a.ab_name < b.ab_name
;


CREATE VIEW IF NOT EXISTS dms_mab_view
AS
SELECT *
FROM dms_single_mab_view
UNION
SELECT *
FROM dms_combo_mab_view
;
