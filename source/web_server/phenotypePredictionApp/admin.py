from django.contrib import admin
import os

# Register your models here.
from .models import Job, PicaModel, Bin


class JobDisplay(admin.ModelAdmin):

    def bal_acc_cutoff(obj):
        return obj.requested_balac
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
            return os.path.basename(str(obj.fileInput))
        return obj.fileInput
    def has_add_permissions(self, request):
        return False
    def has_delete_permissions(self, request, obj=None):
        return False
    no_errors.boolean = True
    link_to_results.allow_tags = True
    actions = None
    list_display = ('job_date', 'user_email', 'user_ip', 'finished_bins', 'total_bins', bal_acc_cutoff, "error_type", uploaded_file_name, link_to_results)
    list_display_links = None
    search_fields = ('user_email', 'user_ip', 'filename')

class ModelDisplay(admin.ModelAdmin):

    def has_add_permissions(self, request):
        return False
    def has_delete_permissions(self, request, obj=None):
        return False
    list_display = ('model_name', 'model_desc', 'model_train_date')
    list_display_links = None
    actions = None

class BinDisplay(admin.ModelAdmin):

    def has_add_permissions(self, request):
        return False
    def has_delete_permissions(self, request, obj=None):
        return False
    def completeness(obj):
        return obj.comple
    def contamination(obj):
        return obj.conta

    list_display = ('bin_name', completeness, contamination, "strainhet", "tax_id", "taxon_name", "taxon_rank")


admin.site.site_header = "PhenDB Admin Pages"
admin.site.site_title = "PhenDB Admin Pages"
admin.site.index_title = "PhenDB Admin Pages"
admin.site.register(Job, JobDisplay)
admin.site.register(PicaModel, ModelDisplay)
admin.site.register(Bin, BinDisplay)
