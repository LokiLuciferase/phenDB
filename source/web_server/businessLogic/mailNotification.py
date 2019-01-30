from django.core.mail import send_mail
import time
import os
import threading
from django.core.mail import *
import subprocess
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
            obj = Job.objects.get(key=self.key)
            if obj.total_bins == obj.finished_bins and obj.total_bins != 0:
                self.__sendMail(obj.user_email, "https://phen.csb.univie.ac.at/phendb/results/" + self.key)
                break
            sleepTime = self.initialSleep * self.runCounter
            if sleepTime > self.maxSleep:
                sleepTime = self.maxSleep
            time.sleep(sleepTime)

    def __sendMail(self, mailAddress, url):
        # TODO: had to set /var/spool/clientmqueue to 777 to allow sending by httpd
        # Fix this soon
        ps = Popen(["/usr/sbin/sendmail", mailAddress], stdin=PIPE, stderr=PIPE)

        message = EmailMessage()
        message.from_email = "donotreply@phen.csb.univie.ac.at"
        message.subject = "PhenDB notification"
        message.body = 'Your PhenDB results are now available under ' + url + '\n \nThis mail was sent automatically. Please do not respond to it.'

        ps.stdin.write(message.message().as_bytes())
        (stdout, stderr) = ps.communicate()

        with open("/apps/phenDB/logs/logmail.txt", "w") as file_log:
            file_log.write("stdout: ")
            file_log.write(stdout)
            file_log.write("stderr: ")
            file_log.write(stderr)
