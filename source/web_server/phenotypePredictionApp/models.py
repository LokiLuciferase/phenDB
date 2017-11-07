from django.db import models
import uuid

# Create your models here.

def upload_function(instance):
    subfolder = instance.key
    filename = instance.filename
    return "documents/" + subfolder + "/" + filename

#not working, TODO: repair
class UploadedFile(models.Model):
    key = models.TextField(default=uuid.uuid4())
    filename = models.TextField()
    file_input = models.FileField(upload_to = upload_function)
    job_date = models.DateTimeField(auto_now_add=True)

class ResultFile(models.Model):
    actualID = models.TextField
    document = models.FileField(upload_to=('resultFiles/' + str(actualID) + '/'))
