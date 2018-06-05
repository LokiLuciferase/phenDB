/*get top @LIM enogs corresponding with YES for model @MODEL_NAME*/
SELECT enog_name, enog_descr, internal_rank from phenotypePredictionApp_enogrank as enogrank
  JOIN phenotypePredictionApp_enog as enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = (SELECT id FROM phenotypePredictionApp_picamodel
                                                  WHERE model_name = 'AEROBE')
ORDER BY internal_rank DESC
LIMIT 10;

/*get top @LIM enogs corresponding with NO for model @MODEL_NAME*/
SELECT enog_name, enog_descr, internal_rank from phenotypePredictionApp_enogrank as enogrank
  JOIN phenotypePredictionApp_enog as enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = (SELECT id FROM phenotypePredictionApp_picamodel
                                                  WHERE model_name = 'AEROBE')
ORDER BY internal_rank ASC
LIMIT 10;