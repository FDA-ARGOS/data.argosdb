CREATE VIEW IF NOT EXISTS rx_dms_single_mab_view
AS
SELECT *
FROM rx_single_mab_view
UNION
SELECT *
FROM dms_single_mab_view
;


CREATE VIEW IF NOT EXISTS rx_dms_combo_mab_view
AS
SELECT *
FROM rx_combo_mab_view
UNION
SELECT *
FROM dms_combo_mab_view
;

CREATE VIEW IF NOT EXISTS rx_dms_mab_view
AS
SELECT *
FROM rx_dms_single_mab_view
UNION
SELECT *
FROM rx_dms_combo_mab_view
;
