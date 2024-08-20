# -*- coding: utf-8 -*-
"""R-Kit model fields

Copyright (c) 2012-13 Christopher S. Case and David H. Bronke
Licensed under the MIT license; see the LICENSE file for details.

"""
from collections import MutableSet

from riak.mapreduce import RiakLink
from riak.riak_object import RiakObject

import rkit
from .validators import validate_number, validate_string, validate_bool, validate_choices


class ModelField(object):
    validators = list()
    index_suffix = "_bin"

    def __init__(self, name=None, default=None, index=False, null=False, choices=None, searchable=False,
            required=False, primary=False):
        self.default = default
        self.index = index
        self.nullable = null
        self.choices = choices
        self.searchable = searchable
        self.required = required
        self.primary = primary

        self.name = name

    def init(self, cls):
        pass

    def __get__(self, instance, type):
        return instance._values.get(self.name, self.default)

    def __set__(self, instance, value):
        instance._values[self.name] = value

        if self.index:
            self.addSecondaryIndex(instance, value)

        if self.primary:
            instance.pk = unicode(value).encode('utf-8')

    def validate(self, instance):
        value = instance._values.get(self.name, self.default)

        # Handle choices support
        if value is not None or self.required:
            validate_choices(self.choices, value)

        # Handle validator support
        for validator in self.validators:
            validator(value)

    def addSecondaryIndex(self, instance, value):
        if isinstance(value, (list, tuple)):
            idx = [unicode(item) for item in value]
        else:
            idx = unicode(value)
        instance._indexes[self.name] = idx

    def renderValue(self, instance):
        return instance._values.get(self.name, self.default)


class NumberField(ModelField):
    validators = [validate_number]
    index_suffix = "_int"

    def __set__(self, instance, value):
        value = float(value)

        super(NumberField, self).__set__(instance, value)


class StringField(ModelField):
    validators = [validate_string]

    def __set__(self, instance, value):
        value = unicode(value)

        super(StringField, self).__set__(instance, value)


class BooleanField(ModelField):
    validators = [validate_bool]

    def __set__(self, instance, value):
        value = bool(value)

        super(BooleanField, self).__set__(instance, value)


class ListField(ModelField):
    def __set__(self, instance, value):
        if not isinstance(value, (list, tuple)):
            raise ValueError("Not a valid list or tuple!")

        value = list(value)

        super(ListField, self).__set__(instance, value)


class DictField(ModelField):
    #TODO: Implement some basic checking here.
    pass


class ObjectCollection(MutableSet):
    """Represents a set of objects. Objects may be added as RiakLinks, RiakObjects, "bucket/key" strings, (bucket, key)
    tuples or lists, and Model instances.

    """
    def __init__(self, objects=set()):
        self._bucket_keys = set()  # A set of (bucket, key) tuples for all objects in the collection.
        self._objects = dict()  # A dictionary from (bucket, key) to Model instance for all instances we already have.

        if objects:
            self += objects

    def __contains__(self, value):
        bucketKey = self._getBucketKey(value)

        return bucketKey in self._bucket_keys

    def __iter__(self):
        for bucket, key in self._bucket_keys:
            yield self._itemFor(bucket, key)

    def __len__(self):
        return len(self._bucket_keys)

    def add(self, value):
        from .base import Model

        bucketKey = self._getBucketKey(value)

        self._bucket_keys.add(bucketKey)

        if isinstance(value, Model):
            #TODO: Subscribe to key/bucket update events on value so we can update self._bucket_keys and self._objects!
            self._objects[bucketKey] = value

    def discard(self, value):
        from .base import Model

        if isinstance(value, Model):
            self._objects.discard(value)

        bucketKey = self._getBucketKey(value)

        self._bucket_keys.discard(bucketKey)

        try:
            del self._objects[bucketKey]
        except KeyError:
            pass

    def get(self):
        if len(self) > 1:
            raise ValueError("LinkField returned more than one value.")

        # Return the first (only) item we have.
        for bucket, key in self._bucket_keys:
            return self._itemFor(bucket, key)

    def bucket_keys(self):
        """Iterate over all contained objects, returning a (bucket, key) tuple for each.

        """
        return iter(self._bucket_keys)

    def links(self):
        """Iterate over all contained objects, returning a RiakLink for each.

        """
        for bucket, key in self._bucket_keys:
            yield RiakLink(bucket, key)

    def _getBucketKey(self, value):
        from .base import Model

        if isinstance(value, Model):
            return value.Meta.bucket, value.pk

        elif isinstance(value, (list, tuple)):
            return value[0], value[1]

        elif isinstance(value, basestring):
            try:
                bucket, key = value.split('/')
                return bucket, key
            except ValueError:
                pass

        elif isinstance(value, (RiakLink, RiakObject)):
            return value.get_bucket(), value.get_key()

        # If we got here, something is very, very wrong.
        raise ValueError("{!r} is not a valid object or link specification!".format(value))

    def _itemFor(self, bucket, key):
        try:
            return self._objects[bucket, key]
        except KeyError:
            pass

        from .base import Model

        cls = Model.bucketModels.get(bucket, None)
        if cls is None:
            #TODO: Return a generic model class instead of a string
            #TODO: Will the base to_python work?
            return "{}/{}".format(bucket, key)

        else:
            # Query the object's data from Riak, and then instantiate a new Model instance for it.
            data = rkit.connection.bucket(bucket).get(key)
            return cls.to_python(data)

    def __iadd__(self, itemOrItems):
        """The '+=' operator, adding either an individual item or all items from an iterable. (treated as iterable by
        default)

        """
        try:
            self |= itemOrItems
        except TypeError:
            self.add(itemOrItems)

    def __isub__(self, itemOrItems):
        """Override the '-=' operator to take either an individual item, or an iterable of items.

        """
        try:
            self.remove(itemOrItems)
        except KeyError:
            self ^= itemOrItems


class LinkCollection(ModelField):
    one_to_one = False

    def __init__(self, related_model=None, related_name="", *args, **kwargs):

        self.related_name = related_name
        self.related_model = related_model

        super(LinkCollection, self).__init__(*args, **kwargs)

    def init(self, cls):
        if self.related_model is not None:

            if isinstance(self.related_model, basestring):
                self.related_model = __import__(cls.__module__).__dict__.get(self.related_model)

            related_name = self.related_name
            if related_name in [None, ""]:
                related_name = cls.Meta.bucket

            if not hasattr(cls, related_name):
                field = LinkCollection(name=related_name)
                self.related_model._fields[related_name] = field
                setattr(self.related_model, related_name, field)

    def __get_ObjectCollection(self, instance):
        links = instance._links.get(self.name)

        if links is None:
            links = ObjectCollection()
            instance._links[self.name] = links

        return links

    def __get__(self, instance, type):
        return self.__get_ObjectCollection(instance)

    def __set__(self, instance, value):
        links = self.__get_ObjectCollection(instance)
        links.clear()
        links += value

    def renderValue(self, instance):
        return None


class LinkField(LinkCollection):
    def __init__(self, one_to_one=False, *args, **kwargs):

        self.one_to_one = one_to_one

        super(LinkField, self).__init__(*args, **kwargs)

    def init(self, cls):
        if self.related_model is not None:

            if isinstance(self.related_model, basestring):
                self.related_model = __import__(cls.__module__).__dict__.get(self.related_model)

            related_name = self.related_name
            if related_name in (None, ""):
                related_name = cls.Meta.bucket

            if not hasattr(cls, related_name):
                field = LinkCollection(name=related_name)

                # If we're a one_to_one field, we then make this a link field, instead of a link collection
                if self.one_to_one:
                    field = LinkField(name=related_name)

                self.related_model._fields[related_name] = field
                setattr(self.related_model, related_name, field)

    def __get__(self, instance, type):
        links = super(LinkField, self).__get__(instance, type)

        return links.get()


class SecondaryIndexField(ModelField):
    def __get__(self, instance, type):
        return instance._indexes[self.name]

    def __set__(self, instance, value):
        self.addSecondaryIndex(instance, value)

    def renderValue(self, instance):
        return None
