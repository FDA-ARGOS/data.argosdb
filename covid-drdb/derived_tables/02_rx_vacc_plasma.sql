INSERT INTO rx_vacc_plasma
  SELECT
    SbjP.ref_name,
    SbjP.rx_name,
    SbjP.subject_name,
    SbjInf.infected_var_name,
    (
      SELECT STRING_AGG(DISTINCT vaccine_name, ' + ')
      FROM subject_vaccines SV
      WHERE
        SV.ref_name=SbjP.ref_name AND
        SV.subject_name=SbjP.subject_name AND
        SV.vaccination_date<SbjP.collection_date
    ) AS vaccine_name,
    SbjP.location,
    GREATEST(ROUND((SbjP.collection_date - SbjVacc.vaccination_date) / 30.), 1) AS timing,
    SbjVacc.dosage,
    SbjP.collection_date,
    SbjP.cumulative_group
  FROM
    subject_plasma SbjP
  JOIN subject_vaccines SbjVacc ON
    SbjVacc.ref_name = SbjP.ref_name AND
    SbjVacc.subject_name = SbjP.subject_name AND
    SbjP.collection_date > SbjVacc.vaccination_date
  LEFT JOIN subject_infections SbjInf ON
    SbjInf.ref_name = SbjP.ref_name AND
    SbjInf.subject_name = SbjP.subject_name AND
    infection_date <= SbjP.collection_date AND
    NOT EXISTS (
      SELECT 1 FROM subject_infections SbjInfNext
      WHERE
        SbjInfNext.ref_name = SbjP.ref_name AND
        SbjInfNext.subject_name = SbjP.subject_name AND
        SbjInfNext.infection_date <= SbjP.collection_date AND
        SbjInf.infection_date < SbjInfNext.infection_date
    )
  WHERE
    NOT EXISTS (
      SELECT 1 FROM subject_vaccines SbjVaccNext
      WHERE
        SbjVaccNext.ref_name = SbjP.ref_name AND
        SbjVaccNext.subject_name = SbjP.subject_name AND
        SbjP.collection_date > SbjVaccNext.vaccination_date AND
        SbjVaccNext.vaccination_date > SbjVacc.vaccination_date
    );

UPDATE rx_vacc_plasma SET timing=NULL WHERE timing=0;
