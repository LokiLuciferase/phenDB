
SET @MODEL_ID = 47;  # sulfate reducer

/* get all PhenDB models including their optimal mean balanced accuracy */
SELECT model_name as "Model Name",
  model_desc AS "Model Description",
  model_train_date AS "Model Train Date",
  ROUND(mean_balanced_accuracy, 2) as "Maximal Accuracy"
  FROM phenotypePredictionApp_picamodel as picamodel
  JOIN phenotypePredictionApp_picamodelaccuracy AS pma
    ON pma.model_id = picamodel.id
      WHERE comple = 1
      AND conta = 0;

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

/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SELECT assembly_id AS "Genome Accession",
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