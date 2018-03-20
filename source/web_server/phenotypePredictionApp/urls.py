from django.conf.urls import url, handler403, handler404, handler500

from . import views

app_name = 'phenotypePredictionApp'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'sendinput$', views.sendinput, name='sendinput'),
    url(r'results/[\S]{36}/$', views.getResults, name="getResults"),
    url(r'results/[\S]{36}/download/$', views.fileDownload, name="fileDownload")
]

handler403 = 'phenotypePredictionApp.views.permissionDenied'
handler404 = 'phenotypePredictionApp.views.pageNotFound'
handler500 = 'phenotypePredictionApp.views.serverError'