CREATE VIEW IF NOT EXISTS rx_conv_plasma_wt_infect_view
AS
SELECT
    *
FROM
    rx_conv_plasma
WHERE
    infected_var_name IN (
        SELECT var_name FROM isolate_wildtype_view
    )
    OR
    infected_var_name IS NULL
;

CREATE VIEW IF NOT EXISTS rx_conv_plasma_infect_var_view
AS
SELECT
    rx.*,
    var.var_name,
    var.as_wildtype
FROM
    rx_conv_plasma rx
LEFT JOIN
    variants var
ON
    rx.infected_var_name = var.var_name
;
