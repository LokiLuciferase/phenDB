#code from https://djangosnippets.org/snippets/201/

from django.core.serializers import serialize
from django.db.models.query import QuerySet
import json
from django.template import Library

register = Library()

def jsonify(object):
    testDir = {"A" : "1"}
    if isinstance(object, QuerySet):
        #for singleObj in object:
            #singleObj.bin.bin_name
        return serialize('json', testDir)
    return serialize('json', testDir)

register.filter('jsonify', jsonify)