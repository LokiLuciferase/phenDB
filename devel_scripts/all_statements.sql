/* get all current PhenDB models including their optimal mean balanced accuracy */
SELECT picamodel.id as 'model_id', model_name as 'model_name',
  model_desc AS 'model_description',
  max(model_train_date) as 'model_train_date',
  ROUND(mean_balanced_accuracy, 2) as 'max_accuracy'
  FROM phenotypePredictionApp_picamodel as picamodel
  JOIN phenotypePredictionApp_picamodelaccuracy AS pma
    ON pma.model_id = picamodel.id
      WHERE comple = 1
      AND conta = 0
  GROUP BY model_name;


/*get training data (accession ID an their associated verdict) for model with id @MODEL_ID*/
SET @MODEL_ID = 47;  # sulfate reducer
SELECT assembly_id,
  taxon_name AS 'scientific_name',
  CASE WHEN verdict = 0 THEN 'NO' ELSE 'YES' END 'trait_presence'
FROM phenotypePredictionApp_picamodeltrainingdata AS picamodeltrainingdata
  JOIN phenotypePredictionApp_taxon AS taxon ON picamodeltrainingdata.tax_id = taxon.tax_id
WHERE picamodeltrainingdata.model_id = @MODEL_ID;


/*get 100 most influential enogs for model with id @MODEL_ID*/
SET @MODEL_ID = 47;  # sulfate reducer
SELECT enog_name AS 'enog_name',
  enog_descr AS 'enog_description',
  internal_rank AS 'rank_in_model',
  ROUND(score, 6) AS 'weight_in_model'
FROM phenotypePredictionApp_enogrank AS enogrank
  JOIN phenotypePredictionApp_enog AS enog ON enogrank.enog_id = enog.id
WHERE enogrank.model_id = @MODEL_ID
ORDER BY internal_rank ASC
LIMIT 100;


/* view a list of refseq genomes calculated for a given model, linked in via model ID from overview.
   Respects balanced accuracy, pica confidence and hierarchical filtering cutoffs*/
SET @MODEL_ID = 47;  # sulfate reducer
SET @ACC_CO =  0.8; # MBA cutoff
SET @CONF_CO = 0.8; # PICA confidence cutoff
SELECT bin.assembly_id as 'assembly_id',
  bin.tax_id as 'taxon_id',
  taxon.taxon_name as 'scientific_name',
  CASE WHEN pr.verdict = 0 THEN 'NO' ELSE 'YES' END 'prediction',
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


/*get a list of refseq pre-calculated genomes with a field which can link to the results via the bin_id*/
SELECT bin.id as 'bin_id',
  bin.assembly_id as 'assembly_id',
  bin.tax_id as 'taxon_id',
  taxon.taxon_name as 'scientific_name',
  'View' as 'link_to_predictions' from phenotypePredictionApp_bininjob as bininjob
  JOIN phenotypePredictionApp_bin as bin on bininjob.bin_id = bin.id
  JOIN phenotypePredictionApp_taxon as taxon on bin.tax_id = taxon.tax_id
  WHERE bin.assembly_id IS NOT NULL
  AND bininjob.job_id = (select job.id from phenotypePredictionApp_job as job
                           WHERE job.key='PHENDB_PRECALC');


/* view a list of traits for phenDB precalculated genome, linked in via bin ID from last query.
   Respects balanced accuracy, pica confidence and hierarchical filtering cutoffs*/
SET @BIN_ID = 14; # one of precalculated refseq genomes
SET @ACC_CO = 0.8; # MBA cutoff
SET @CONF_CO = 0.8; # PICA confidence cutoff
SELECT model.id as 'model_id',#
       model.model_name AS 'trait_name',
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
  GROUP BY model.model_name
  ORDER BY model.model_train_date;
