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

    def get_absolute_url(self):
        return "results/%s/" % self.key


class Bin(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['md5sum'])
        ]
        verbose_name = 'PhenDB bin'
        verbose_name_plural = 'PhenDB bins'

    bin_name = models.TextField()
    md5sum = models.CharField(unique=True, max_length=32)
    tax_id = models.TextField(null=True, blank=True)
    taxon_name = models.TextField(default="")
    taxon_rank = models.TextField(default="")
    comple = models.FloatField()
    conta = models.FloatField()
    strainhet = models.FloatField()
    assembly_id = models.TextField(null=True, blank=True)  # TODO: migrate this!

    def __str__(self):
        return str(self.bin_name) + " " + str(self.md5sum) + " " + str(self.tax_id) + " " + str(self.taxon_name) + " " + str(self.taxon_rank) + " " + str(self.comple) + " " + str(self.conta) + " " + str(self.strainhet)
        #return "File name: {fn}".format(fn=self.bin_name)

    def natural_key(self):
        return (self.bin_name)


class BinInJob(models.Model):

    class Meta:
        unique_together = ('bin', 'job')  # composite primary key
        indexes = [
            models.Index(fields=['bin', 'job']),
            models.Index(fields=['bin',]),
            models.Index(fields=['job',])
       ]

    bin = models.ForeignKey(Bin, on_delete=models.CASCADE)  # delete bins if we delete the association row
    job = models.ForeignKey(Job)

    def __str__(self):
        return "Bin: {bin} in Job {job}".format(bin=self.bin, job=self.job)


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


class PicaModel(models.Model):

    class Meta:
        unique_together = ('model_name', 'model_train_date')  # composite primary key
        indexes = [
            models.Index(fields=['model_name', 'model_train_date'])
        ]
        verbose_name = 'PICA Model'
        verbose_name_plural = 'PICA Models'

    # New fields must be added at end due to web UI
    model_name = models.CharField(max_length=64)
    type = models.CharField(max_length=2)
    model_desc = models.TextField()
    model_train_date = models.DateTimeField(auto_now=True)

    def __str__(self):
        return "Name: {mid}\t Description: {md}\tDate of Training: {mtd} \t " \
               "Type= {type}".format(mid=self.model_name,md=self.model_desc, mtd=str(self.model_train_date),
                                         type=self.type)

    def natural_key(self):
        return (self.model_name)


class EnogRank(models.Model):
#todo: change PK to model and enog alone
    class Meta:
        unique_together = ('model', 'enog')  # composite primary key
        indexes = [
            models.Index(fields=['model', 'enog']),
            models.Index(fields=['model',]),
            models.Index(fields=['enog',]),
            models.Index(fields=['internal_rank',])
        ]

    model = models.ForeignKey(PicaModel)
    enog = models.ForeignKey(Enog)
    internal_rank = models.FloatField()

    def __str__(self):
        return "Enog ID: {eid}\tModel ID: {mid}\tRank {ir}\t(Model trained on {mtd})".format(eid=self.enog_id,
                                                                                           mid=self.model.model_name,
                                                                                           ir=str(self.internal_rank),
                                                                                           mtd=str(self.model.model_train_date))


class PicaModelAccuracy(models.Model):
    class Meta:
        unique_together = ('model', 'comple', 'conta')  # composite primary key
        indexes = [
            models.Index(fields=['model', 'comple', 'conta'])
        ]

    model = models.ForeignKey(PicaModel)
    comple = models.FloatField()
    conta = models.FloatField()
    mean_balanced_accuracy  = models.FloatField()
    mean_fp_rate = models.FloatField()
    mean_fn_rate = models.FloatField()


    def __str__(self):
        return "Name: {mid}\t type= {type} \t Completeness: {comple} \t " \
               "Contamination {conta} \t bal. acc {balac}".format(mid=self.model.model_name,
                                                                      type=self.model.type, comple=self.comple,
                                                                      conta=self.conta,
                                                                      balac=self.mean_balanced_accuracy
                                                                              )

class PicaModelTrainingData(models.Model):
    class Meta:
        unique_together = ('model', 'tax_id', 'assembly_id')
        indexes = [
            models.Index(fields=['model', 'tax_id'])
        ]

    model = models.ForeignKey(PicaModel)
    tax_id = models.CharField(max_length=10)
    assembly_id = models.CharField(max_length=20)
    verdict = models.NullBooleanField()

    def __str__(self):
        return "Taxid {ti} used in model {mod} (counted as {yn})".format(ti=self.taxid,
                                                                         mod=self.model.model_name,
                                                                         yn=self.verdict)


class HmmerResult(models.Model):

    class Meta:
        unique_together = ('bin', 'enog')
        indexes = [
            models.Index(fields=['enog', 'bin']),
            models.Index(fields=['bin'])
        ]

    enog = models.ForeignKey(Enog)
    bin = models.ForeignKey(Bin)

    def __str__(self):
        return "Enog {eid} contained in the bin {mds} of Job {Jb}".format(eid=self.enog_id,
                                                              mds=self.bin_id, jb=self.bin.key)


class PicaResult(models.Model):

    class Meta:
        unique_together = ('bin', 'model')
        indexes = [
            models.Index(fields=['bin', 'model']),
            models.Index(fields=['bin'])
        ]

    bin = models.ForeignKey(Bin)
    model = models.ForeignKey(PicaModel)
    verdict = models.NullBooleanField()
    pica_pval = models.FloatField()
    accuracy = models.FloatField()

    def __str__(self):
        return "Bin {mds}: {v} ({acc} accuracy, {pval} pica_p_value) for model {mid} trained on {mtd}.".format(mds=self.bin.bin_name,
                                                                                          v=str(self.verdict),
                                                                                          acc=str(self.accuracy),
                                                                                        pval=str(self.pica_pval),
                                                                                          mid=self.model.model_name,
                                                                                          mtd=str(self.model.model_train_date))


class Taxon(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['tax_id', ]),
        ]

    tax_id = models.CharField(max_length=10)
    taxon_name = models.TextField()

    def __str__(self):
        return "{tn} - {ti}".format(tn=self.taxon_name, ti=self.tax_id)