CREATE VIEW IF NOT EXISTS rx_single_mab_view
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

CREATE VIEW IF NOT EXISTS rx_combo_mab_view
AS
SELECT
    ref_name,
    rx_name,
    GROUP_CONCAT(ab_name, '/'),
    '' AS synonyms,
    'mixed' AS availability,
    '' AS pdb_id,
    '' AS target,
    '' AS class,
    '' AS epitope,
    '' AS institute,
    '' AS origin
FROM rx_antibodies
GROUP BY ref_name, rx_name
HAVING count(1) > 1
ORDER BY ab_name;


CREATE VIEW IF NOT EXISTS rx_mab_view
AS
SELECT *
FROM rx_single_mab_view
UNION
SELECT *
FROM rx_combo_mab_view
;
