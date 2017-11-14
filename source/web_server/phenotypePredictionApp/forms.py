from django import forms
from .models import UploadedFile


class FileForm(forms.ModelForm):
    class Meta:
        model = UploadedFile
        exclude = ()
    def __init__(self, *args, **kwargs):
        super(FileForm, self).__init__(*args, **kwargs)
        self.fields['key'].initial = self.data['key']
        self.fields['filename'].initial = self.data['filename']