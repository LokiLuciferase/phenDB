# -*- coding: utf-8 -*-
# Generated by Django 1.11.7 on 2017-11-19 21:41
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion
import phenotypePredictionApp.models
import uuid


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='bin',
            fields=[
                ('bin_id', models.TextField()),
                ('file_name', models.TextField()),
                ('genome_path', models.TextField()),
                ('md5sum', models.TextField(primary_key=True, serialize=False)),
                ('errors', models.TextField()),
                ('comple', models.FloatField()),
                ('conta', models.FloatField()),
            ],
        ),
        migrations.CreateModel(
            name='enog',
            fields=[
                ('enog_id', models.TextField(primary_key=True, serialize=False)),
                ('enog_descr', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='job',
            fields=[
                ('job_id', models.TextField(default=uuid.UUID('8761b815-e7f6-4eda-88da-9615676a3a6a'), primary_key=True, serialize=False)),
                ('user_ip', models.TextField()),
                ('user_email', models.TextField()),
                ('job_date', models.DateTimeField(auto_now=True)),
                ('folder_path', models.TextField()),
                ('output_tgz', models.FileField(upload_to=phenotypePredictionApp.models.upload_function_upload)),
                ('job_status', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='model',
            fields=[
                ('model_name', models.TextField()),
                ('version_nr', models.IntegerField()),
                ('mname_vnr', models.TextField(primary_key=True, serialize=False)),
                ('is_newest', models.BooleanField()),
                ('model_desc', models.TextField()),
                ('model_train_date', models.DateField(auto_now=True)),
            ],
        ),
        migrations.CreateModel(
            name='model_enog_ranks',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('internal_rank', models.FloatField()),
                ('enog', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.enog')),
                ('model', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.model')),
            ],
        ),
        migrations.CreateModel(
            name='result_enog',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('bin', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.bin')),
                ('enog', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.enog')),
            ],
        ),
        migrations.CreateModel(
            name='result_model',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('verdict', models.BooleanField()),
                ('accuracy', models.FloatField()),
                ('bin', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.bin')),
                ('model', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.model')),
            ],
        ),
        migrations.CreateModel(
            name='ResultFile',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('actualID', models.TextField()),
                ('document', models.FileField(upload_to=phenotypePredictionApp.models.upload_function_results)),
            ],
        ),
        migrations.CreateModel(
            name='UploadedFile',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('key', models.TextField(default=uuid.UUID('41d14c3c-808f-48c0-9a21-5c703701183b'))),
                ('filename', models.TextField()),
                ('fileInput', models.FileField(upload_to=phenotypePredictionApp.models.upload_function_upload)),
            ],
        ),
        migrations.AddIndex(
            model_name='job',
            index=models.Index(fields=['job_id'], name='phenotypePr_job_id_8dfaf8_idx'),
        ),
        migrations.AddIndex(
            model_name='enog',
            index=models.Index(fields=['enog_id'], name='phenotypePr_enog_id_11be84_idx'),
        ),
        migrations.AddField(
            model_name='bin',
            name='job',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='phenotypePredictionApp.job'),
        ),
        migrations.AddIndex(
            model_name='result_model',
            index=models.Index(fields=['bin', 'model'], name='phenotypePr_bin_id_422581_idx'),
        ),
        migrations.AlterUniqueTogether(
            name='result_model',
            unique_together=set([('bin', 'model')]),
        ),
        migrations.AddIndex(
            model_name='result_enog',
            index=models.Index(fields=['enog', 'bin'], name='phenotypePr_enog_id_0f1063_idx'),
        ),
        migrations.AlterUniqueTogether(
            name='result_enog',
            unique_together=set([('bin', 'enog')]),
        ),
        migrations.AddIndex(
            model_name='model_enog_ranks',
            index=models.Index(fields=['model', 'enog'], name='phenotypePr_model_i_85ef6b_idx'),
        ),
        migrations.AlterUniqueTogether(
            name='model_enog_ranks',
            unique_together=set([('model', 'enog', 'internal_rank')]),
        ),
        migrations.AddIndex(
            model_name='bin',
            index=models.Index(fields=['md5sum'], name='phenotypePr_md5sum_0ceaa4_idx'),
        ),
    ]
