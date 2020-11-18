from django.db import models
import uuid
from django.contrib.contenttypes.models import ContentType


# ------------ functions for file upload -----------------------
def upload_function_upload(instance, filename):
    print('upload called')
    subfolder = instance.key
    filename = instance.filename
    return "uploads/" + subfolder + "/" + filename


def upload_function_results(instance, givenname):
    print(givenname)
    subfolder = instance.key
    filename = "results/" + subfolder + "_results/" + givenname
    return filename

# ----------------Models--------------------------------------
# Contains information on the job submitted and file uploaded by the user,
# as well as the output file to be accessed.
class Job(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['key'])
        ]
        verbose_name = 'PhenDB job'
        verbose_name_plural = 'PhenDB jobs'

    key = models.CharField(default=uuid.uuid4(), unique=True, max_length=36)
    filename = models.TextField()
    fileInput = models.FileField(upload_to=upload_function_upload)
    fileOutput = models.FileField(upload_to=upload_function_results)
    user_ip = models.TextField(null=True, blank=True)
    user_email = models.TextField(null=True, blank=True)
    job_date = models.DateTimeField(auto_now=True)
    errors = models.NullBooleanField(null=True, blank=True)
    error_type = models.TextField(default="")
    finished_bins = models.IntegerField(default='0')
    total_bins = models.IntegerField(default='0')
    requested_balac = models.FloatField(default='0.5')
    requested_conf = models.FloatField(default='0.5')
    disable_cutoffs = models.BooleanField(default=False)
    get_explanations = models.BooleanField(default=False)


    def get_absolute_url(self):
        return "results/%s/" % self.key


# Contains information on uploaded genomes.
# Uniquely identified by md5sum hash of genome file.
class Bin(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['md5sum'])
        ]
        verbose_name = 'PhenDB bin'
        verbose_name_plural = 'PhenDB bins'

    md5sum = models.CharField(unique=True, max_length=32)
    tax_id = models.TextField(null=True, blank=True)
    assembly_id = models.TextField(null=True, blank=True)
    comple = models.FloatField()
    conta = models.FloatField()
    strainhet = models.FloatField()


    def __str__(self):
        return " ".join([str(x) for x in (self.md5sum, self.tax_id, self.comple,
                                          self.conta, self.strainhet)])


# Correlates Jobs with Bins; each bin has an alias
# which is the name of the file when it was submitted with the given job
class BinInJob(models.Model):

    class Meta:
        unique_together = ('bin', 'job')  # composite primary key
        indexes = [
            models.Index(fields=['bin', 'job']),
            models.Index(fields=['bin',]),
            models.Index(fields=['job',])
        ]

    bin = models.ForeignKey(Bin, on_delete=models.CASCADE)
    job = models.ForeignKey(Job, on_delete=models.CASCADE)
    bin_alias = models.TextField(null=True, blank=True)

    def __str__(self):
        return "Bin: {bin} in Job {job}".format(bin=self.bin, job=self.job)


# Contains information on EggNOG Enogs (clusters of orthologous groups) which are the features used by PICA
class Enog(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['enog_name'])
        ]

    enog_name = models.CharField(unique=True, max_length=32)
    enog_descr = models.TextField()

    def __str__(self):
        return "ID: {eid}\tDescription: {ed}".format(eid=self.enog_name,
                                                     ed=self.enog_descr)


# Contains information on PICA models used for predicting traits,
# uniquely identified by name of the model and date of training.
class PicaModel(models.Model):

    class Meta:
        unique_together = ('model_name', 'model_train_date')
        indexes = [
            models.Index(fields=['model_name', 'model_train_date'])
        ]
        verbose_name = 'PICA Model'
        verbose_name_plural = 'PICA Models'

    # New fields must be added at end due to web UI
    model_name = models.CharField(max_length=64)
    alg_type = models.CharField(max_length=64)
    feature_type = models.CharField(max_length=64)
    model_desc = models.TextField()
    model_train_date = models.DateTimeField(auto_now=True)

    def __str__(self):
        return "Name: {mid}\t Description: {md}\tDate of Training: {mtd} \t " \
               "Type= {type}".format(mid=self.model_name, md=self.model_desc,
                                     mtd=str(self.model_train_date), type=self.type)


# Contains the ranking of PICA model features (meaning,
# the importance of any enog in the model for the decision of that model).
class EnogRank(models.Model):

    class Meta:
        unique_together = ('model', 'enog')
        indexes = [
            models.Index(fields=['model', 'enog']),
            models.Index(fields=['model',]),
            models.Index(fields=['enog',]),
            models.Index(fields=['internal_rank',])
        ]

    model = models.ForeignKey(PicaModel, on_delete=models.CASCADE)
    enog = models.ForeignKey(Enog, on_delete=models.CASCADE)
    internal_rank = models.IntegerField()
    score = models.FloatField()
    pred_class = models.BooleanField()

    def __str__(self):
        return "Enog ID: {eid}\tModel ID:" \
               " {mid}\tRank {ir}\t(Model trained on {mtd})".format(eid=self.enog_id,
                                                                    mid=self.model.model_name,
                                                                    ir=str(self.internal_rank),
                                                                    mtd=str(self.model.model_train_date))


# Contains precomputed estimated accuracy values of any PICA model
# at a certain completeness and contamination level of the bin in question.
class PicaModelAccuracy(models.Model):
    class Meta:
        unique_together = ('model', 'comple', 'conta')
        indexes = [
            models.Index(fields=['model', 'comple', 'conta'])
        ]

    model = models.ForeignKey(PicaModel, on_delete=models.CASCADE)
    comple = models.FloatField()
    conta = models.FloatField()
    mean_balanced_accuracy  = models.FloatField()
    mean_fp_rate = models.FloatField(null=True)
    mean_fn_rate = models.FloatField(null=True)

    def __str__(self):
        return "Name: {mid}\t type= {type} \t Completeness: {comple} \t " \
               "Contamination {conta} \t bal. acc {balac}".format(mid=self.model.model_name,
                                                                  type=self.model.type, comple=self.comple,
                                                                  conta=self.conta,
                                                                  balac=self.mean_balanced_accuracy)


# Contains training instances for each PICA model,
# incl. taxon id, assembly id and the class of the training instance (YES or NO)
class PicaModelTrainingData(models.Model):
    class Meta:
        unique_together = ('model', 'tax_id', 'assembly_id')
        indexes = [
            models.Index(fields=['model', 'tax_id'])
        ]

    model = models.ForeignKey(PicaModel, on_delete=models.CASCADE)
    tax_id = models.CharField(max_length=10)
    assembly_id = models.CharField(max_length=64)
    verdict = models.NullBooleanField()

    def __str__(self):
        return "Taxid {ti} used in model {mod} (counted as {yn})".format(ti=self.tax_id,
                                                                         mod=self.model.model_name,
                                                                         yn=self.verdict)


# Contains Enogs present in a given bin, as found by HMMER.
class HmmerResult(models.Model):

    class Meta:
        unique_together = ('bin', 'enog')
        indexes = [
            models.Index(fields=['enog', 'bin']),
            models.Index(fields=['bin'])
        ]

    enog = models.ForeignKey(Enog, on_delete=models.CASCADE)
    bin = models.ForeignKey(Bin, on_delete=models.CASCADE)

    def __str__(self):
        return "Enog {eid} contained in the bin {mds}.".format(eid=self.enog_id,
                                                               mds=self.bin_id)


# Contains predictions of trait presence or absence for a given bin and model, including confidence scores
class PicaResult(models.Model):

    class Meta:
        unique_together = ('bin', 'model')
        indexes = [
            models.Index(fields=['bin', 'model']),
            models.Index(fields=['bin'])
        ]

    bin = models.ForeignKey(Bin, on_delete=models.CASCADE)
    model = models.ForeignKey(PicaModel, on_delete=models.CASCADE)
    verdict = models.NullBooleanField()
    pica_pval = models.FloatField()
    accuracy = models.FloatField()
    nc_masked = models.BooleanField(default=False)

    def __str__(self):
        return "Bin {mds}: {v} ({acc} accuracy, " \
               "{pval} pica_p_value) for model {mid} trained on {mtd}.".format(mds=self.bin.md5sum,
                                                                               v=str(self.verdict),
                                                                               acc=str(self.accuracy),
                                                                               pval=str(self.pica_pval),
                                                                               mid=self.model.model_name,
                                                                               mtd=str(self.model.model_train_date))


class PicaResultExplanation(models.Model):

    class Meta:
        unique_together = ('bin', 'model', 'enog')
        indexes = [
            models.Index(fields=['bin', 'model', 'enog']),
            models.Index(fields=['bin', 'model']),
            models.Index(fields=['bin'])
        ]

    bin  = models.ForeignKey(Bin, on_delete=models.CASCADE)
    model = models.ForeignKey(PicaModel, on_delete=models.CASCADE)
    enog = models.ForeignKey(Enog, on_delete=models.CASCADE)
    enog_is_present = models.BooleanField()
    delta_shap = models.FloatField()

    def __str__(self):
        return "Bin {mds}, Model {mo}: presence={enog_is_present} of Enog {e}" \
               " ==> delta SHAP of {ds}.".format(
            mds=str(self.bin.md5sum),
            mo=str(self.model.model_name),
            enog_is_present=str(self.enog_is_present),
            e=str(self.enog.enog_name),
            ds=str(self.delta_shap)
        )


# Contains a stripped down version of the NCBI Taxonomy Names table,
# which only contains Scientific Names and their correlated Taxon IDs.
class Taxon(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['tax_id', ]),
        ]

    tax_id = models.CharField(unique=True, max_length=10)
    taxon_rank = models.TextField(default="")
    taxon_name = models.TextField()

    def __str__(self):
        return "{tn} - {ti}".format(tn=self.taxon_name, ti=self.tax_id)
