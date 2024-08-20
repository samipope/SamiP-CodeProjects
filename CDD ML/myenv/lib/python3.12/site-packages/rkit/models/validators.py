# -*- coding: utf-8 -*-
"""R-Kit validators

Copyright (c) 2012-13 Christopher S. Case and David H. Bronke
Licensed under the MIT license; see the LICENSE file for details.

"""
from .exceptions import ValidationException


def validate_choices(choices, value):
    if choices is not None:
        try:
            if not isinstance(value, (list, tuple)):
                value = [value]

            for val in value:
                if not val in choices:
                    raise ValidationException("Value not a valid choice!")
        except ValidationException as ex:
            raise ex
        except:
            raise ValidationException("Error evaluating choices.")


def validate_number(value):
    try:
        float(value)
    except:
        raise ValidationException("Not a valid number!")


def validate_string(value):
    try:
        unicode(value)
    except:
        raise ValidationException("Not a valid string!")


def validate_bool(value):
    try:
        bool(value)
    except:
        raise ValidationException("Not a valid bool!")
