from django.db import models
import uuid
from django.contrib.contenttypes.models import ContentType


#------------ functions for file upload -----------------------
def upload_function_upload(instance, filename):
    print('upload called')
    subfolder = instance.key
    filename = instance.filename
    return "documents/" + subfolder + "/" + filename

def upload_function_results(instance, givenname):
    print(givenname)
    subfolder = instance.key
    filename = "resultFiles/" + subfolder + "_results/" + givenname
    return filename


#----------------Models--------------------------------------

#TODO: rename later (in the end)
class UploadedFile(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['key'])
        ]

    key = models.CharField(default=uuid.uuid4(), unique=True, max_length=36)
    filename = models.TextField()
    fileInput = models.FileField(upload_to=upload_function_upload)
    fileOutput = models.FileField(upload_to=upload_function_results)
    user_ip = models.TextField(null=True, blank=True)
    user_email = models.TextField(null=True, blank=True)
    job_date = models.DateTimeField(auto_now=True)
    folder_path = models.TextField()
    errors = models.NullBooleanField(null=True, blank=True)
    job_status = models.TextField(default='0')
    def get_absolute_url(self):
        return "results/%s/" % self.key


class bin(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['md5sum'])
        ]

    bin_name = models.TextField()
    md5sum = models.CharField(unique=True, max_length=32)
    comple = models.FloatField()
    conta = models.FloatField()

    def __str__(self):
        return "File name: {fn}".format(fn=self.bin_name)


class bins_in_UploadedFile(models.Model):

    class Meta:
        unique_together = ('bin', 'UploadedFile') # composite primary key
        indexes = [
            models.Index(fields=['bin', 'UploadedFile'])
       ]

    bin = models.ForeignKey(bin)
    UploadedFile=models.ForeignKey(UploadedFile)

    def __str__(self):
        return "Bin: {bin} in Job {job}".format(bin=self.bin, job=self.UploadedFile)


class enog(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['enog_name'])
        ]

    enog_name = models.CharField(unique=True, max_length=32)
    enog_descr = models.TextField()

    def __str__(self):
        return "ID: {eid}\tDescription: {ed}".format(eid=self.enog_name,
                                                     ed=self.enog_descr)


class model(models.Model):

    class Meta:
        unique_together = ('model_name', 'model_train_date') # composite primary key
        indexes = [
            models.Index(fields=['model_name', 'model_train_date'])
       ]

    model_name = models.CharField(max_length=64)
    type = models.CharField(max_length=2)
    model_desc = models.TextField()
    model_train_date = models.DateTimeField(auto_now=True)

    def __str__(self):
        return "Name: {mid}\t Description: {md}\tDate of Training: {mtd} \t " \
               "Type= {type}".format(mid=self.model_name,md=self.model_desc, mtd=str(self.model_train_date),
                                         type=self.type)


class model_enog_ranks(models.Model):
#todo: change PK to model and enog alone
    class Meta:
        unique_together = ('model', 'enog')  # composite primary key
        indexes = [
            models.Index(fields=['model', 'enog'])
        ]

    model = models.ForeignKey(model)
    enog = models.ForeignKey(enog)
    internal_rank = models.FloatField()

    def __str__(self):
        return "Enog ID: {eid}\tModel ID: {mid}\tRank {ir}\t(Model trained on {mtd})".format(eid=self.enog_id,
                                                                                           mid=self.model.model_name,
                                                                                           ir=str(self.internal_rank),
                                                                                           mtd=str(self.model.model_train_date))
class model_accuracies(models.Model):
    class Meta:
        unique_together = ('model', 'comple', 'conta')  # composite primary key
        indexes = [
            models.Index(fields=['model', 'comple', 'conta'])
        ]

    model = models.ForeignKey(model)
    comple = models.FloatField()
    conta = models.FloatField()
    mean_balanced_accuracy  = models.FloatField()
    mean_fp_rate = models.FloatField()
    mean_fn_rate = models.FloatField()


    def __str__(self):
        return "Name: {mid}\t type= {type} \t Completness: {comple} \t " \
               "Contamination {conta} \t bal. acc {balac}".format(mid=self.model.model_name,
                                                                      type=self.model.type, comple=self.comple,
                                                                      conta=self.conta,
                                                                      balac=self.mean_balanced_accuracy
                                                                              )

class result_enog(models.Model):

    class Meta:
        unique_together = ('bin', 'enog')
        indexes = [
            models.Index(fields=['enog', 'bin'])
        ]

    enog = models.ForeignKey(enog)
    bin = models.ForeignKey(bin)

    def __str__(self):
        return "Enog {eid} contained in the bin {mds} of Job {Jb}".format(eid=self.enog_id,
                                                              mds=self.bin_id, jb=self.UploadedFile.key)


class result_model(models.Model):

    class Meta:
        unique_together = ('bin', 'model')
        indexes = [
            models.Index(fields=['bin', 'model'])
        ]

    bin = models.ForeignKey(bin)
    model = models.ForeignKey(model)
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
