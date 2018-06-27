

/* get all current PhenDB models including their optimal mean balanced accuracy */
SELECT model_name as "Model Name",
  model_desc AS "Model Description",
  max(model_train_date) AS "Model Train Date",
  ROUND(mean_balanced_accuracy, 2) as "Maximal Accuracy"
  FROM phenotypePredictionApp_picamodel as picamodel
  JOIN phenotypePredictionApp_picamodelaccuracy AS pma
    ON pma.model_id = picamodel.id
      WHERE comple = 1
      AND conta = 0
  GROUP BY model_name;

SET @MODEL_ID = 47;  # sulfate reducer
/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SELECT assembly_id AS "Genome Accession",
  taxon_name AS "Scientific Name",
  CASE WHEN verdict = 0 THEN "NO" ELSE "YES" END "Trait Presence"
FROM phenotypePredictionApp_picamodeltrainingdata AS picamodeltrainingdata
  JOIN phenotypePredictionApp_taxon AS taxon ON picamodeltrainingdata.tax_id = taxon.tax_id
WHERE picamodeltrainingdata.model_id = @MODEL_ID;

SET @MODEL_ID = 47;  # sulfate reducer
/*get 100 most influential enogs corresponding with NO for model with id @MODEL_ID*/
SELECT enog_name AS "Enog Name",
  enog_descr AS "Enog Description",
  internal_rank AS "Rank in Model",
  ROUND(score, 6) AS "Weight in Model"
FROM phenotypePredictionApp_enogrank AS enogrank
  JOIN phenotypePredictionApp_enog AS enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
  AND pred_class = 0
ORDER BY internal_rank ASC
LIMIT 100;

SET @MODEL_ID = 47;  # sulfate reducer
/*get 100 most influential enogs corresponding with YES for model with id @MODEL_ID*/
SELECT enog_name AS "Enog Name",
  enog_descr AS "Enog Description",
  internal_rank AS "Rank in Model",
  ROUND(score, 6) AS "Weight in Model"
FROM phenotypePredictionApp_enogrank AS enogrank
  JOIN phenotypePredictionApp_enog AS enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
  AND pred_class = 1
ORDER BY internal_rank ASC
LIMIT 100;

SET @MODEL_ID = 47;  # sulfate reducer
SET @ACC_CO =  0.8; # MBA cutoff
SET @CONF_CO = 0.8; # PICA confidence cutoff
/* view a list of refseq genomes calculated for a given model, linked in via model ID from overview.
   Respects balanced accuracy, pica confidence and hierarchical filtering cutoffs*/
SELECT bin.assembly_id as 'assembly_id',
  taxon.taxon_name as 'scientific_name',
  CASE WHEN pr.verdict = 0 THEN "NO" ELSE "YES" END 'prediction',
  ROUND(pr.pica_pval, 2) AS 'prediction_confidence',
  ROUND(pr.accuracy, 2) AS 'balanced_accuracy' from phenotypePredictionApp_bin AS bin
  JOIN phenotypePredictionApp_picaresult as pr ON bin.id = pr.bin_id
  JOIN phenotypePredictionApp_taxon as taxon on bin.tax_id = taxon.tax_id
  WHERE pr.model_id = @MODEL_ID
  AND pr.accuracy > @ACC_CO
  AND pr.pica_pval > @CONF_CO
  AND pr.nc_masked = 0
  AND bin.assembly_id IS NOT NULL
  AND bin.id IN (select bin_id from phenotypePredictionApp_bininjob as bij
                           WHERE bij.job_id = (select job.id from phenotypePredictionApp_job as job
                                               WHERE job.key='PHENDB_PRECALC'));


/*get a list of refseq pre-calculated genomes with a field which can link to the results*/
SELECT bin.assembly_id as 'assembly_id',
  taxon.taxon_name as 'scientific_name',
  'View' as 'link_to_predictions' from phenotypePredictionApp_bininjob as bininjob
  JOIN phenotypePredictionApp_bin as bin on bininjob.bin_id = bin.id
  JOIN phenotypePredictionApp_taxon as taxon on bin.tax_id = taxon.tax_id
  WHERE bin.assembly_id IS NOT NULL
  AND bininjob.job_id = (select job.id from phenotypePredictionApp_job as job
                           WHERE job.key='PHENDB_PRECALC');


SET @BIN_ID = 14; # one of precalculated refseq genomes
SET @ACC_CO = 0.8; # MBA cutoff
SET @CONF_CO = 0.8; # PICA confidence cutoff
/* view a list of traits for phenDB precalculated genome, linked in via bin ID from last query.
   Respects balanced accuracy, pica confidence and hierarchical filtering cutoffs*/
SELECT model.model_name AS 'trait_name',
       CASE
       WHEN pr.nc_masked = 0 THEN (CASE
                                   WHEN pr.pica_pval < @CONF_CO OR pr.accuracy < @CONF_CO
                                   THEN 'ND'
                                   ELSE (CASE
                                         WHEN pr.verdict = 0
                                         THEN 'NO'
                                         ELSE 'YES'
                                         END)
                                   END)
       ELSE 'NC'
       END 'prediction',
       ROUND(pr.pica_pval, 2) AS 'prediction_confidence',
       ROUND(pr.accuracy, 2) AS 'balanced_accuracy',
       model.model_desc as 'model_description' FROM phenotypePredictionApp_picaresult AS pr
  JOIN phenotypePredictionApp_picamodel AS model ON pr.model_id = model.id
  WHERE pr.bin_id = @BIN_ID
  ORDER BY model.model_train_date;
