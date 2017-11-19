from django.db import models
import uuid
from django.contrib.contenttypes.models import ContentType

def upload_function(instance):
    subfolder = instance.job_id
    filename = instance.file_name
    return "documents/" + subfolder + "/" + filename

# LL: the following two tables should be removed after functionality of web interface has been ported
# to accomodate the database strcture defined further down


class UploadedFile(models.Model):
    key = models.TextField(default=uuid.uuid4())
    filename = models.TextField()
    file_input = models.FileField(upload_to = upload_function)
    job_date = models.DateTimeField(auto_now_add=True)


class ResultFile(models.Model):
    actualID = models.TextField
    document = models.FileField(upload_to=('resultFiles/' + str(actualID) + '/'))


class job(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['job_name'])
        ]

    job_name = models.TextField(default=uuid.uuid4(), primary_key=True)
    user_ip = models.TextField()
    user_email = models.TextField()
    job_date = models.DateTimeField(auto_now=True)
    folder_path = models.TextField()
    output_tgz = models.FileField(upload_to=upload_function)
    job_status = models.TextField()

    def __str__(self):
        return "User IP: {ip}\nRecieved: {jd}".format(ip=self.user_ip, jd=self.job_date)


class bin(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['md5sum'])
        ]

    bin_name = models.TextField()
    job = models.ForeignKey(job)
    md5sum = models.TextField(primary_key=True)
    errors = models.TextField()
    comple = models.FloatField()
    conta = models.FloatField()

    def __str__(self):
        return "File name: {fn}\nErrors: {err}".format(fn=self.bin_name,
                                                       err=self.errors if self.errors else "")


class enog(models.Model):

    class Meta:
        indexes = [
            models.Index(fields=['enog_name'])
        ]

    enog_name = models.TextField(primary_key=True)
    enog_descr = models.TextField()
    enog_category = models.TextField()

    def __str__(self):
        return "ID: {eid}\tDescription: {ed}\tCategor(y/ies): {ca}".format(eid=self.enog_name,
                                                     ed=self.enog_descr, ca=self.enog_category)


class model(models.Model):

    class Meta:
        unique_together = ('model_name', 'version_nr') # composite primary key
        indexes = [
            models.Index(fields=['model_name', 'is_newest'])
       ]

    model_name = models.TextField()
    version_nr = models.IntegerField()
    is_newest = models.BooleanField()
    model_desc = models.TextField()
    model_train_date = models.DateField(auto_now=True)

    def __str__(self):
        return "Name: {mid}\tVersion Nr: {vnr}\t Description: {md}\tDate of Training: {mtd} \t is_newest= {new}".format(mid=self.model_name,
                                                                                md=self.model_desc,
                                                                                mtd=str(self.model_train_date), vnr=self.version_nr, new=self.is_newest)


class model_enog_ranks(models.Model):

    class Meta:
        unique_together = ('model', 'enog', 'internal_rank')  # composite primary key
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


class result_enog(models.Model):

    class Meta:
        unique_together = ('bin', 'enog')
        indexes = [
            models.Index(fields=['enog', 'bin'])
        ]

    enog = models.ForeignKey(enog)
    bin = models.ForeignKey(bin)

    def __str__(self):
        return "Enog {eid} contained in the bin {mds}".format(eid=self.enog_id,
                                                              mds=self.bin_id)


class result_model(models.Model):

    class Meta:
        unique_together = ('bin', 'model')
        indexes = [
            models.Index(fields=['bin', 'model'])
        ]

    bin = models.ForeignKey(bin)
    model = models.ForeignKey(model)
    verdict = models.BooleanField()
    accuracy = models.FloatField()

    def __str__(self):
        return "Bin {mds}: {v} ({acc} accuracy) for model {mid} trained on {mtd}.".format(mds=self.bin.bin_name,
                                                                                          v=str(self.verdict),
                                                                                          acc=str(self.accuracy),
                                                                                          mid=self.model.model_name,
                                                                                          mtd=str(self.model.model_train_date))
