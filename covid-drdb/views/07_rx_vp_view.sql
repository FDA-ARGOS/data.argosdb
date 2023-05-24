CREATE VIEW IF NOT EXISTS rx_vacc_plasma_infect_var_view
AS
SELECT
    rx.*,
    var.var_name,
    var.as_wildtype
FROM
    rx_vacc_plasma rx
LEFT JOIN
    variants var
ON
    rx.infected_var_name = var.var_name
;


-- CREATE VIEW IF NOT EXISTS rx_vacc_plasma_full_dosage_view
-- AS
-- SELECT
--     rx.*
-- FROM
--     rx_vacc_plasma rx,
--     vaccines vac
-- WHERE
--     rx.vaccine_name = vac.vaccine_name
--     AND
--     rx.dosage = vac.st_shot
-- ;
