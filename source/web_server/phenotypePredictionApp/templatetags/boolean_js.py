from django.template import Library

register = Library()

@register.filter(is_safe=True)
def boolean_js(object):
    if(object):
        return "true"
    else:
        return "false"