/*get training data (accession ID an their associated verdict) for model @MODEL_NAME*/
SELECT assembly_id, tax_id, verdict FROM phenotypePredictionApp_picamodeltrainingdata as picamodeltrainingdata
  JOIN phenotypePredictionApp_picamodel as picamodel ON picamodeltrainingdata.model_id = picamodel.id
WHERE picamodel.model_name = 'AEROBE';