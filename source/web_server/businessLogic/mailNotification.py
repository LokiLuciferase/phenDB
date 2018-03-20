from django.core.mail import send_mail
import time
import threading
from django.core.mail import *
import subprocess
from phenotypePrediction.settings import GlobalVariables
from phenotypePredictionApp.models import *
from subprocess import Popen, PIPE


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

        #message =  'To:' + mailAddress + '\n Subject: phenDB notification \n From: donotreply@phen.csb.univie.ac.at \n Your phenDB results are now available under phen.csb.univie.ac.at' + url + '\n \n This mail was sent automatically.Please do not respond to it.'

        ps = Popen(["/usr/sbin/sendmail"] + mailAddress, stdin=PIPE, stderr=PIPE)

        message = EmailMessage()
        message.to = mailAddress
        message.subject = "phenDB notification"
        message.body = 'Your phenDB results are now available under phen.csb.univie.ac.at' + url + '\n \n This mail was sent automatically.Please do not respond to it.'

        ps.stdin.write(message.message().as_bytes())
        (stdout, stderr) = ps.communicate()

        print("mailNotification called")
        print(stdout)
        print(stderr)

        file_log = open("/apps/phenDB/logs/logmail.txt", "w")
        file_log.write("stdout:")
        file_log.write(stdout)
        file_log.write("stderr:")
        file_log.write(stderr)

