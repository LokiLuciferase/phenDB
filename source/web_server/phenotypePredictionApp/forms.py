from django import forms
from .models import UploadedFile


class FileForm(forms.ModelForm):
    class Meta:
        model = UploadedFile
        fields = ['filename', 'key', 'fileInput', 'user_ip', 'user_email']
    def __init__(self, *args, **kwargs):
        super(FileForm, self).__init__(*args, **kwargs)
        self.fields['key'].initial = self.data['key']
        self.fields['filename'].initial = self.data['filename']
        self.fields['fileInput'].initial = self.data['fileInput']
        self.fields['user_ip'].initial = self.data['user_ip']
        self.fields['user_email'].initial = self.data['user_email']
        for key, value in self.data.items():
            print(str(key) + ' '  + str(value))