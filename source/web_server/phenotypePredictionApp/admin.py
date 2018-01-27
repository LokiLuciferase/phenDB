from django.contrib import admin

# Register your models here.

from .models import UploadedFile

class JobDisplay(admin.ModelAdmin):
    list_display = ('job_date', 'user_email', 'user_ip', 'finished_bins', 'total_bins', 'errors', 'fileInput', 'fileOutput')

admin.site.register(UploadedFile, JobDisplay)