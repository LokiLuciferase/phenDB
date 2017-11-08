# -*- coding: utf-8 -*-
# Generated by Django 1.11.7 on 2017-11-08 15:24
from __future__ import unicode_literals

from django.db import migrations, models
import uuid


class Migration(migrations.Migration):

    dependencies = [
        ('phenotypePredictionApp', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='job',
            name='job_id',
            field=models.TextField(default=uuid.UUID('cab169e0-dcc8-4868-a666-eca873bd2f4c'), primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='uploadedfile',
            name='key',
            field=models.TextField(default=uuid.UUID('b8f5f807-7d78-43b0-8063-4d568c65feec')),
        ),
    ]
