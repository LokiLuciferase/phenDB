

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
/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SELECT assembly_id AS "Genome Accession",
  taxon_name AS "Scientific Name",
  CASE WHEN verdict = 0 THEN "NO" ELSE "YES" END "Trait Presence"
FROM phenotypePredictionApp_picamodeltrainingdata AS picamodeltrainingdata
  JOIN phenotypePredictionApp_taxon AS taxon ON picamodeltrainingdata.tax_id = taxon.tax_id
WHERE picamodeltrainingdata.model_id = @MODEL_ID;

# SET @MODEL_ID = 47;  # sulfate reducer
# /*get mean balanced accuracy, mean fn and fp rate at the desired completeness/contamination for model with id @MODEL_ID*/
# SELECT comple AS "Completeness",
#   conta AS "Contamination",
#   mean_balanced_accuracy AS "Mean Balanced Accuracy"
# FROM phenotypePredictionApp_picamodelaccuracy AS picamodelaccuracy
# WHERE picamodelaccuracy.model_id = @MODEL_ID
#   AND comple >= 0.9
#   AND conta <= 0.1
#   ORDER BY comple DESC
# LIMIT 10;

/*get a list of refseq pre-calculated genomes with a field which can link to the results*/
SELECT bin.taxon_name as "Scientific Name",
  bin.assembly_id as "Assembly Accession",
  'View' as "Link to Trait Predictions" from phenotypePredictionApp_bininjob as bininjob
  JOIN phenotypePredictionApp_bin as bin on bininjob.bin_id = bin.id
  WHERE bininjob.job_id = (select job.id from phenotypePredictionApp_job as job
                           WHERE job.key='PHENDB_PRECALC');


SET @BIN_ID = 14; # one of precalculated refseq genomes
/* view a list of traits for phenDB precalculated genome, linked in via bin ID from last query */
SELECT model.model_name AS "Trait Name",
       CASE WHEN pr.verdict = 0 THEN "NO" ELSE "YES" END "Prediction",
       ROUND(pr.pica_pval, 2) AS "Model Confidence",
       ROUND(pr.accuracy, 2) AS "Balanced Accuracy",
       model.model_desc as "Model Description" FROM phenotypePredictionApp_picaresult AS pr
  JOIN phenotypePredictionApp_picamodel AS model ON pr.model_id = model.id
  WHERE pr.bin_id = @BIN_ID
  ORDER BY model.model_train_date;


SET @MODEL_ID = 47;  # sulfate reducer
/* view a list of refseq genomes calculated for a given model, linked in via model ID from overview */
SELECT bin.assembly_id as "Assembly Accession",
  bin.taxon_name as "Scientific Name",
  CASE WHEN pr.verdict = 0 THEN "NO" ELSE "YES" END "Prediction",
  ROUND(pr.pica_pval, 2) AS "Model Confidence",
  ROUND(pr.accuracy, 2) AS "Balanced Accuracy" from phenotypePredictionApp_bin AS bin
  JOIN phenotypePredictionApp_picaresult as pr ON bin.id = pr.bin_id
  WHERE pr.model_id = @MODEL_ID
  #AND bin.assembly_id IS NOT NULL ## activate this once taxid updating after precalc works
  AND bin.id IN (select bin_id from phenotypePredictionApp_bininjob as bij
                           WHERE bij.job_id = (select job.id from phenotypePredictionApp_job as job
                                               WHERE job.key='PHENDB_PRECALC'));