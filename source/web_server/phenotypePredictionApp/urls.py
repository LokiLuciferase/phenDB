from django.conf.urls import url
from django.conf.urls import handler403, handler404, handler500

from . import views

app_name = 'phenotypePredictionApp'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'sendinput$', views.sendinput, name='sendinput'),
    url(r'results/[\S]{36}/$', views.getResults, name="getResults"),
    url(r'results/[\S]{36}/download/$', views.fileDownload, name="fileDownload")
]

handler403 = 'views.permissionDenied'
handler404 = 'views.pageNotFound'
handler500 = 'views.serverError'