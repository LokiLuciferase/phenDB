from django.contrib import admin
import os

# Register your models here.

from .models import UploadedFile


class JobDisplay(admin.ModelAdmin):

    def no_errors(obj):
        if obj.errors == None:
            return None
        else:
            return not obj.errors

    def link_to_results(obj):
        if obj.fileOutput:
            return "phen.csb.univie.ac.at/phendb/results/{k}".format(k=obj.key)
        return obj.fileOutput

    def uploaded_file_name(obj):
        if obj.fileInput:
            return os.path.basename(obj.fileInput)
        return obj.fileInput

    no_errors.boolean = True
    link_to_results.allow_tags = True

    list_display = ('job_date', 'user_email', 'user_ip', 'finished_bins', 'total_bins', no_errors, uploaded_file_name, link_to_results)
    list_display_links = None

admin.site.register(UploadedFile, JobDisplay)