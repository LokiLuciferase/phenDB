/*get mean balanced accuracy, mean fn and fp rate at the desired completeness/contamination for model @MODEL_NAME*/
SELECT comple, conta, mean_balanced_accuracy, mean_fn_rate, mean_fp_rate FROM phenotypePredictionApp_picamodelaccuracy AS picamodelaccuracy
  JOIN phenotypePredictionApp_picamodel as picamodel ON picamodelaccuracy.model_id = picamodel.id
WHERE picamodel.model_name = 'AEROBE'
  AND comple >= 0.9
  AND conta <= 0.1
  ORDER BY comple DESC
LIMIT 10;