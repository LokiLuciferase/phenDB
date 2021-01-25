import time
import threading
from django.core.mail import EmailMessage
from phenotypePredictionApp.models import Job
from subprocess import Popen, PIPE


class MailNotification(threading.Thread):
    @staticmethod
    def send_mail_now(mail_address: str, content: str, subject: str = 'PhenDB Dev Notice') -> int:
        # TODO: had to set /var/spool/clientmqueue to 777 to allow sending by httpd
        ps = Popen(["/usr/sbin/sendmail", mail_address], stdin=PIPE, stderr=PIPE)
        message = EmailMessage()
        message.from_email = "donotreply@phen.csb.univie.ac.at"
        message.subject = subject
        message.body = content
        ps.stdin.write(message.message().as_bytes())
        (stdout, stderr) = ps.communicate()
        return ps.returncode

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
                body = f"Your PhenDB results are now available under" \
                       f" https://phen.csb.univie.ac.at/phendb/results/{self.key}" \
                       f"\n\nThis mail was sent automatically. Please do not respond to it."
                MailNotification.send_mail_now(
                    mail_address=obj.user_email,
                    content=body,
                    subject='PhenDB notification'
                )
                break
            sleepTime = self.initialSleep * self.runCounter
            if sleepTime > self.maxSleep:
                sleepTime = self.maxSleep
            time.sleep(sleepTime)
