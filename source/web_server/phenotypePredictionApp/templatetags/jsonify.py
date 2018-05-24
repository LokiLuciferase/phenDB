#code from https://djangosnippets.org/snippets/201/

from django.core.serializers import serialize
from django.db.models.query import QuerySet
import json
from django.template import Library

register = Library()

@register.filter(is_safe=True)
def jsonify(object):
    print("filter type object: " + type(object))
    print("filter size object: " + len(object))
    if isinstance(object, QuerySet):
        #for singleObj in object:
            #singleObj.bin.bin_name
        return json.dumps(object)
    return json.dumps(object)