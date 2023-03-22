CREATE VIEW IF NOT EXISTS susc_results_50_omicron_view
AS
SELECT
    *
FROM
    susc_results_50_view susc,
    isolate_omicron_view omi
WHERE
    susc.control_iso_name = omi.iso_name
;


CREATE VIEW IF NOT EXISTS susc_results_50_ba2_view
AS
SELECT
    *
FROM
    susc_results_50_view susc,
    isolate_omicron_ba2_view omi
WHERE
    susc.control_iso_name = omi.iso_name
;
