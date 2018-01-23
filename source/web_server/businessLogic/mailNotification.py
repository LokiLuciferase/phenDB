from django.core.mail import send_mail
import time
import threading
from phenotypePrediction.settings import GlobalVariables
from phenotypePredictionApp.models import *


class MailNotification(threading.Thread):

    def __init__(self, key):
        threading.Thread.__init__(self)
        self.key = key
        self.runCounter = 0
        self.initialSleep = 10  # initial time for thread to sleep in seconds, will be increased dynamically (linear to runCounter), max:20min
        self.maxSleep = 1200

    def run(self):
        while True:
            self.runCounter += 1
            obj = UploadedFile.objects.get(key=self.key)
            if obj.total_bins == obj.finished_bins and obj.total_bins != 0:
                self.__sendMail(obj.user_email, GlobalVariables.WEBSERVER_URL + "/phendb/results/" + self.key)
                break
            sleepTime = self.initialSleep * self.runCounter
            if sleepTime > self.maxSleep:
                sleepTime = self.maxSleep
            time.sleep(sleepTime)

    def __sendMail(self, mailAddress, url):
        send_mail(
            'phenDB notification',
            'Your phenDB results are now available under ' + url + '\n \n This mail was sent automatically. Please do not respond to it.',
            'webapptest@gmx.de',
            [mailAddress],
            fail_silently=True)
