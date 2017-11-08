# -*- coding: utf-8 -*-
# Generated by Django 1.11.7 on 2017-11-08 15:25
from __future__ import unicode_literals

from django.db import migrations, models
import uuid


class Migration(migrations.Migration):

    dependencies = [
        ('phenotypePredictionApp', '0002_auto_20171108_1624'),
    ]

    operations = [
        migrations.AlterField(
            model_name='job',
            name='job_id',
            field=models.TextField(default=uuid.UUID('8b6244d9-eb67-4217-9f3b-7f46eabca7f7'), primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='uploadedfile',
            name='key',
            field=models.TextField(default=uuid.UUID('77c244ac-a52b-4e94-b93c-23fd8b6066a4')),
        ),
    ]
