from django import forms
from .models import job


class FileForm(forms.ModelForm):
    class Meta:
        model = job
        fields = ['filename', 'job_name', 'upload_path']
    def __init__(self, *args, **kwargs):
        super(FileForm, self).__init__(*args, **kwargs)
        self.fields['job_name'].initial = self.data['job_name']
        self.fields['filename'].initial = self.data['filename']
        self.fields['upload_path'].initial = self.data['upload_path']
        for key, value in self.data.items():
            print(str(key) + ' '  + str(value))