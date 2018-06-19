
SET @MODEL_ID = 47;  # sulfate reducer

/* get all PhenDB models */
SELECT model_name as "Model Name",
  model_desc AS "Model Description",
  model_train_date AS "Model Train Date"
FROM phenotypePredictionApp_picamodel as picamodel;

/*get 10 most influential enogs corresponding with YES for model with id @MODEL_ID*/
SELECT enog_name AS "Enog Name",
  enog_descr AS "Enog Description",
  internal_rank AS "Rank in Model",
  CASE WHEN pred_class = 0 THEN "NO" ELSE "YES" END "Trait Class Indicated"
FROM phenotypePredictionApp_enogrank AS enogrank
  JOIN phenotypePredictionApp_enog AS enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
  AND pred_class = 0
ORDER BY internal_rank ASC
LIMIT 10;

/*get 10 most influential enogs corresponding with YES for model with id @MODEL_ID*/
SELECT enog_name AS "Enog Name",
  enog_descr AS "Enog Description",
  internal_rank AS "Rank in Model",
  CASE WHEN pred_class = 0 THEN "NO" ELSE "YES" END "Trait Class"
FROM phenotypePredictionApp_enogrank AS enogrank
  JOIN phenotypePredictionApp_enog AS enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
  AND pred_class = 1
ORDER BY internal_rank ASC
LIMIT 10;

/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SELECT assembly_id AS "Assembly ID",
  picamodeltrainingdata.tax_id AS "Taxon ID",
  taxon_name AS "Scientific Name",
  CASE WHEN verdict = 0 THEN "NO" ELSE "YES" END "Trait Presence"
FROM phenotypePredictionApp_picamodeltrainingdata AS picamodeltrainingdata
  JOIN phenotypePredictionApp_taxon AS taxon ON picamodeltrainingdata.tax_id = taxon.tax_id
WHERE picamodeltrainingdata.model_id = @MODEL_ID;

/*get mean balanced accuracy, mean fn and fp rate at the desired completeness/contamination for model with id @MODEL_ID*/
SELECT comple AS "Completeness",
  conta AS "Contamination",
  mean_balanced_accuracy AS "Mean Balanced Accuracy"
FROM phenotypePredictionApp_picamodelaccuracy AS picamodelaccuracy
WHERE picamodelaccuracy.model_id = @MODEL_ID
  AND comple >= 0.9
  AND conta <= 0.1
  ORDER BY comple DESC
LIMIT 10;