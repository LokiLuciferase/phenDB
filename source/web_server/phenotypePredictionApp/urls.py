from django.conf.urls import url

from . import views

app_name = 'phenotypePredictionApp'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'sendinput$', views.sendinput, name='sendinput'),
    url(r'results/[\S]{36}/$', views.getResults, name="getResults"),
    url(r'results/PHENDB_PRECALC/$', views.getResults, name="getResults"),
    url(r'results/[\S]{36}/download/$', views.fileDownload, name="fileDownload")
]
