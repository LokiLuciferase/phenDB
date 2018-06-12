
SET @MODEL_ID = 48;  # motility

/*get 10 most influential enogs corresponding with NO for model with id @MODEL_ID*/
SELECT enog_name, enog_descr, internal_rank from phenotypePredictionApp_enogrank as enogrank
  JOIN phenotypePredictionApp_enog as enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
ORDER BY internal_rank DESC
LIMIT 10;

/*get 10 most influential enogs corresponding with YES for model with id @MODEL_ID*/
SELECT enog_name, enog_descr, internal_rank from phenotypePredictionApp_enogrank as enogrank
  JOIN phenotypePredictionApp_enog as enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
ORDER BY internal_rank ASC
LIMIT 10;

/*get mean balanced accuracy, mean fn and fp rate at the desired completeness/contamination for model with id @MODEL_ID*/
SELECT comple, conta, mean_balanced_accuracy, mean_fn_rate, mean_fp_rate FROM phenotypePredictionApp_picamodelaccuracy AS picamodelaccuracy
WHERE picamodelaccuracy.model_id = @MODEL_ID
  AND comple >= 0.9
  AND conta <= 0.1
  ORDER BY comple DESC
LIMIT 10;

/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SELECT assembly_id, picamodeltrainingdata.tax_id, taxon_name, verdict FROM phenotypePredictionApp_picamodeltrainingdata as picamodeltrainingdata
  JOIN phenotypePredictionApp_taxon as taxon on picamodeltrainingdata.tax_id = taxon.tax_id
WHERE picamodeltrainingdata.model_id = @MODEL_ID;
