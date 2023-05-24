CREATE VIEW IF NOT EXISTS susc_results_view
AS
SELECT *
FROM
    (
        SELECT
            ref_name,
            rx_name,
            rx_group,
            control_iso_name,
            iso_name,
            section,
            fold_cmp,
            fold,
            potency_type,
            control_potency,
            potency,
            potency_unit,
            resistance_level,
            ineffective,
            control_cumulative_count,
            cumulative_count,
            control_assay_name,
            assay_name,
            date_added
        FROM susc_results
        WHERE rx_name IS NOT NULL

        UNION

        SELECT
            susc.ref_name,
            unlink.rx_name,
            susc.rx_group,
            susc.control_iso_name,
            susc.iso_name,
            susc.section,
            susc.fold_cmp,
            susc.fold,
            susc.potency_type,
            susc.control_potency,
            unlink.potency,
            susc.potency_unit,
            susc.resistance_level,
            susc.ineffective,
            susc.control_cumulative_count,
            unlink.cumulative_count,
            susc.control_assay_name,
            unlink.assay_name,
            susc.date_added
        FROM
            susc_results susc,
            unlinked_susc_results unlink
        WHERE
            susc.ref_name = unlink.ref_name
            AND
            susc.rx_group = unlink.rx_group
            AND
            susc.iso_name = unlink.iso_name
            AND
            susc.assay_name = unlink.assay_name
            AND
            susc.potency_type = unlink.potency_type
            AND
            susc.rx_name IS NULL
    )
;

CREATE VIEW IF NOT EXISTS susc_results_50_view
AS
SELECT
    *
FROM
    susc_results_view
WHERE
    potency_type IN ('IC50', 'NT50')
;

CREATE VIEW IF NOT EXISTS susc_results_wt_view
AS
SELECT
    *
FROM
    susc_results_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;

CREATE VIEW IF NOT EXISTS susc_results_50_wt_view
AS
SELECT
    *
FROM
    susc_results_50_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;


CREATE VIEW IF NOT EXISTS susc_results_cp_50_wt_view
AS
SELECT
    *
FROM
    susc_results_50_wt_view susc,
    rx_conv_plasma rx
WHERE
    susc.ref_name = rx.ref_name
    AND
    susc.rx_name = rx.rx_name
;

CREATE VIEW IF NOT EXISTS susc_results_mab_50_wt_view
AS
SELECT
    *
FROM
    susc_results_50_wt_view susc,
    rx_mab_view rx
WHERE
    susc.ref_name = rx.ref_name
    AND
    susc.rx_name = rx.rx_name
;


CREATE VIEW IF NOT EXISTS susc_results_vp_50_wt_view
AS
SELECT
    *
FROM
    susc_results_50_wt_view susc,
    rx_vacc_plasma rx,
    vaccines vac
WHERE
    susc.ref_name = rx.ref_name
    AND
    susc.rx_name = rx.rx_name
    AND
    rx.vaccine_name = vac.vaccine_name
;


CREATE VIEW IF NOT EXISTS susc_results_aggr_view
AS
SELECT
    *
FROM
    susc_results_view a
WHERE
    EXISTS (
        SELECT 1
        FROM susc_results_view b
        WHERE
            a.ref_name = b.ref_name
            AND
            a.rx_name = b.rx_name
            AND
            a.iso_name = b.iso_name
            AND
            a.control_iso_name = b.control_iso_name
            AND
            b.cumulative_count > 1
    )
;


CREATE VIEW IF NOT EXISTS susc_results_aggr_50_view
AS
SELECT
    *
FROM
    susc_results_aggr_view a
WHERE
    potency_type IN ('IC50', 'NT50')
;

CREATE VIEW IF NOT EXISTS susc_results_aggr_wt_view
AS
SELECT
    *
FROM
    susc_results_aggr_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;

CREATE VIEW IF NOT EXISTS susc_results_aggr_50_wt_view
AS
SELECT
    *
FROM
    susc_results_aggr_50_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;

CREATE VIEW IF NOT EXISTS susc_results_indiv_view
AS
SELECT
    *
FROM
    susc_results_view a
WHERE
    NOT EXISTS (
        SELECT 1
        FROM susc_results_view b
        WHERE
            a.ref_name = b.ref_name
            AND
            a.rx_name = b.rx_name
            AND
            a.iso_name = b.iso_name
            AND
            a.control_iso_name = b.control_iso_name
            AND
            b.cumulative_count > 1
    )
;


CREATE VIEW IF NOT EXISTS susc_results_indiv_50_view
AS
SELECT
    *
FROM
    susc_results_indiv_view a
WHERE
    potency_type IN ('IC50', 'NT50')
;

CREATE VIEW IF NOT EXISTS susc_results_indiv_wt_view
AS
SELECT
    *
FROM
    susc_results_indiv_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;

CREATE VIEW IF NOT EXISTS susc_results_indiv_50_wt_view
AS
SELECT
    *
FROM
    susc_results_indiv_50_view susc,
    isolate_wildtype_view wt
WHERE
    susc.control_iso_name = wt.iso_name
;
