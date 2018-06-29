from django.template import Library
import json

register = Library()

@register.filter(is_safe=True)
def jsonify(object):
    return json.dumps(object)