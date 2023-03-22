-- In susc_results, one of fold_cmp/fold or resistance_level should not be NULL
ALTER TABLE susc_results
  ADD CONSTRAINT chk_fold_resistance_level CHECK (
    (
      -- case 1, fold_cmp/fold are presented
      fold_cmp IS NOT NULL AND
      fold IS NOT NULL AND
      resistance_level IS NULL
    ) OR
    (
      -- case 2, resistance_level is presented
      fold_cmp IS NULL AND
      fold IS NULL AND
      resistance_level IS NOT NULL
    )
  );

-- In rx_potency, potency must be greater than 0
CREATE FUNCTION checkPotencyNonzero(
  potency NUMERIC(10,3),
  ref_name VARCHAR,
  rx_name VARCHAR,
  iso_name VARCHAR,
  potency_type VARCHAR,
  lower_limit NUMERIC(10,3),
  upper_limit NUMERIC(10,3)
) RETURNS BOOLEAN
AS
$$
DECLARE
  lim NUMERIC(10,3);
  msg VARCHAR;
BEGIN
  IF potency > 0 THEN
    RETURN TRUE;
  END IF;
  lim := lower_limit;
  msg := FORMAT('(potency_lower_limit for %s)', potency_type);
  IF potency_type LIKE 'IC%' THEN
    lim := upper_limit;
    msg := FORMAT('(potency_upper_limit for %s)', potency_type);
  END IF;
  RAISE EXCEPTION E'rx_potency.potency must be greater than 0 (ref_name=\x1b[1m%\x1b[0m rx_name=\x1b[1m%\x1b[0m iso_name=\x1b[1m%\x1b[0m potency=\x1b[1m%\x1b[0m); consider changing it to \x1b[1m%\x1b[0m %', ref_name, rx_name, iso_name, potency, lim, msg;
END;
$$ LANGUAGE PLPGSQL IMMUTABLE;

ALTER TABLE rx_potency
  ADD CONSTRAINT chk_potency_nonzero CHECK (checkPotencyNonzero(potency, ref_name, rx_name, iso_name, potency_type::text, potency_lower_limit, potency_upper_limit));


ALTER TABLE rx_potency
  ADD CONSTRAINT chk_potency_limit_and_unit CHECK (
    (
      -- case 1
      potency_type IN ('NT50', 'NT80', 'NT90', 'NT100') AND
      potency_lower_limit IS NOT NULL AND
      potency_unit IS NULL
    ) OR
    (
      -- case 2
      potency_type IN ('IC50', 'IC80', 'IC90', 'IC100', 'EC50', 'RC50') AND
      potency_upper_limit IS NOT NULL AND
      potency_unit IN ('ng/ml', 'µM', 'pM', 'nM', 'IU/ml')
    ) OR
    (
      -- case 3
      potency_type IN ('NC20', 'NC') AND
      potency_lower_limit IS NOT NULL AND
      potency_unit = 'percent'
    ) OR
    (
      -- case 4
      potency_type = 'Ki' AND
      potency_unit IN ('ng/ml', 'µM', 'pM', 'nM', 'IU/ml')
    ) OR
    (
      -- case 5
      potency_type = 'Kcat/Km' AND
      potency_lower_limit IS NULL AND
      potency_upper_limit IS NULL AND
      potency_unit IN ('s⁻¹µM⁻¹', 's⁻¹M⁻¹')
    )
  );

-- In subject_vaccines, infection_date of later dose must be later than eariler dose
CREATE FUNCTION checkVaccineDosageOrder(
  rname VARCHAR,
  sname VARCHAR,
  dose INT,
  vacc_date DATE
) RETURNS BOOLEAN
AS $$
  SELECT NOT EXISTS (
    SELECT 1 FROM subject_vaccines
    WHERE rname = ref_name AND sname = subject_name AND (
      (dose > dosage AND vacc_date <= vaccination_date) OR
      (dose < dosage AND vacc_date >= vaccination_date)
    )
  )
$$ LANGUAGE SQL;

ALTER TABLE subject_vaccines
  ADD CONSTRAINT chk_vaccine_dosage_order CHECK (checkVaccineDosageOrder(ref_name, subject_name, dosage, vaccination_date));

-- In subject_treatments, start date must be not greater than end date
ALTER TABLE subject_treatments
  ADD CONSTRAINT chk_sbjrx_time_range CHECK (
    start_date <= end_date
  );

CREATE FUNCTION hasPreviousInfection(rname varchar, sname varchar, event_date date) RETURNS boolean AS $$
  SELECT EXISTS (
    SELECT 1 FROM subject_infections
    WHERE rname = ref_name AND sname = subject_name AND event_date >= infection_date
  )
$$ LANGUAGE SQL;

CREATE FUNCTION hasPreviousVaccine(rname varchar, sname varchar, event_date date) RETURNS boolean AS $$
  SELECT EXISTS (
    SELECT 1 FROM subject_vaccines
    WHERE rname = ref_name AND sname = subject_name AND event_date >= vaccination_date
  )
$$ LANGUAGE SQL;

ALTER TABLE subject_plasma
  ADD CONSTRAINT chk_previous_infection_or_vaccine CHECK (
    hasPreviousInfection(ref_name, subject_name, collection_date) OR
    hasPreviousVaccine(ref_name, subject_name, collection_date)
  );

ALTER TABLE subject_isolates
  ADD CONSTRAINT chk_previous_infection CHECK (
    hasPreviousInfection(ref_name, subject_name, collection_date)
  );

ALTER TABLE isolates
  ADD CONSTRAINT no_private_var_name CHECK (
    var_name NOT LIKE '\_%'
  );

ALTER TABLE subject_infections
  ADD CONSTRAINT no_private_infected_var_name CHECK (
    infected_var_name NOT LIKE '\_%'
  );

ALTER TABLE dms_escape_results
  ADD CONSTRAINT max_escape_score_is_one CHECK (
    escape_score <= 1
  )
