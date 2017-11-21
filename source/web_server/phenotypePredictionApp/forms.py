from django import forms
from .models import job


class FileForm(forms.ModelForm):
    class Meta:
        model = job
        fields = ['filename', 'job_name', 'upload_path', 'user_ip', 'user_email', 'job_date', 'job_status']
    def __init__(self, *args, **kwargs):
        super(FileForm, self).__init__(*args, **kwargs)
        self.fields['job_name'].initial = self.data['job_name']
        self.fields['filename'].initial = self.data['filename']
        self.fields['upload_path'].initial = self.data['upload_path']
        self.fields['user_ip'].initial = self.data['user_ip']
        self.fields['user_email'].initial = self.data['user_email']
        self.fields['job_status'].initial = self.data['job_status']

        for key, value in self.data.items():
            print(str(key) + ' '  + str(value))