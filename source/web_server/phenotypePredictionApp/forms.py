from django import forms
from .models import Job


class FileForm(forms.ModelForm):
    class Meta:
        model = Job
        fields = ['filename', 'key', 'fileInput', 'user_ip', 'user_email', 'errors', 'requested_balac']

    def __init__(self, *args, **kwargs):
        super(FileForm, self).__init__(*args, **kwargs)
        self.fields['key'].initial = self.data['key']
        self.fields['filename'].initial = self.data['filename']
        self.fields['fileInput'].initial = self.data['fileInput']
        self.fields['user_ip'].initial = self.data['user_ip']
        self.fields['user_email'].initial = self.data['user_email']
        self.fields['errors'].initial = self.data['errors']
        self.fields['requested_balac'].initial = self.data['requested_balac']
        for key, value in self.data.items():
            print(str(key) + ' ' + str(value))
