from django.conf.urls import url

from . import views

app_name = 'phenotypePredictionApp'
urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'sendinput', views.sendinput, name='sendinput')
]
