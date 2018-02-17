from django.contrib import admin

# Register your models here.

from .models import UploadedFile


class JobDisplay(admin.ModelAdmin):

    def no_errors(obj):
        if obj.errors == None:
            return None
        else:
            return not obj.errors

    no_errors.boolean = True

    list_display = ('job_date', 'user_email', 'user_ip', 'finished_bins', 'total_bins', no_errors, 'fileInput', 'fileOutput')

admin.site.register(UploadedFile, JobDisplay)