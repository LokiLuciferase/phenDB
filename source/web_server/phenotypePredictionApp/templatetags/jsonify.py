
from django.core.serializers import serialize
from django.template import Library

register = Library()

@register.filter(is_safe=True)
def jsonify(object):
    return serialize('json', object, use_natural_foreign_keys=False)