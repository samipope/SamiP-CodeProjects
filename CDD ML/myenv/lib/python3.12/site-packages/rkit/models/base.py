# -*- coding: utf-8 -*-
"""R-Kit base models

Copyright (c) 2012-13 Christopher S. Case and David H. Bronke
Licensed under the MIT license; see the LICENSE file for details.

"""
from contextlib import contextmanager

# Riak python client imports
import riak

# R-kit imports
import rkit
from .fields import ModelField, ObjectCollection
from .query import ModelManager

#----------------------------------------------------------------------------------------------------------------------


class ModelMeta(type):
    """A metaclass for our Models. This is in charge of populating several special fields on the model class, as well
    as properly handling inheritance issues, and properly initializing our ModelField instances.

    """
    # A mapping from bucket name to model class
    bucketModels = dict()

    def __new__(mcs, name, bases, dict):
        # Setup the model fields
        fields = {}
        for key, val in dict.iteritems():
            if isinstance(val, ModelField):
                if key in ["pk", "_values", "_links", "_indexes", "_fields"]:
                    raise AttributeError("Cannot override reserved field name '{}'.".format(key))

                if val.name is None:
                    val.name = key
                fields[val.name] = val

        dict['_fields'] = fields

        # Create the class
        newCls = super(ModelMeta, mcs).__new__(mcs, name, bases, dict)

        # Set the object manager
        newCls.objects = ModelManager(newCls)

        # Meta class handling
        userMeta = dict.get('Meta', object())

        class Meta(object):
            bucket = getattr(userMeta, 'bucket', name.lower())

        newCls.Meta = Meta

        mcs.bucketModels[newCls.Meta.bucket] = newCls

        # Finish field initialization
        for field in newCls._fields.values():
            field.init(newCls)

        return newCls


class Count(object):
    def __init__(self, initial=0):
        self.value = initial

    def __enter__(self):
        self.value += 1

    def __exit__(self, exc_type, exc_value, traceback):
        self.value -= 1

    def __call__(self):
        return self.value


class Model(object):
    """The base class for all models. Handles most of the model behavior, and all of the conversion to/from riak.

    """
    __metaclass__ = ModelMeta

    instances = dict()
    _pk = None

    _unsavedLinks = set()

    countSaveCalls = Count()

    def __init__(self):
        """Model class initialization.

        """
        if rkit.connection is None:
            raise ValueError("Cannot create model without valid connection object. Did you forget to call 'connect'?")

        self._values = dict()
        self._links = dict()
        self._indexes = dict()
        self.pk = None

        # Make a local copy
        self._fields = self._fields.copy()

    @property
    def bucket(self):
        """A read-only property making it easier to access the bucket name for this model.

        """
        return self.Meta.bucket

    @property
    def pk(self):
        """The primary key of this model.

        """
        return self._pk

    @pk.setter
    def pk(self, value):
        """Sets the primary key. Also updates the global registry of all models with the new information.

        **Never override this property.**

        """
        if self._pk != value:
            # Remove our old key
            if (self.bucket, self._pk) in Model.instances.keys():
                del Model.instances[(self.bucket, self._pk)]

            # Update the dict.
            Model.instances[(self.bucket, value)] = self

            # Update the primary key
            self._pk = value

    def __getattr__(self, name):
        """Checks the list of links for the requested attribute, and returns that. This just makes working with links
        more convienent.

        """
        try:
            return self._links[name]
        except KeyError:
            return None

    def save(self):
        """Saves the model to Riak. You may override this function, however, you will either want to still call this
        base method, or you will want to perform the same steps as this:

        * Perform a full_clean.
        * Get this model as a Riak object
        * Save the riak object

        """
        with self.countSaveCalls:
            self.full_clean()

            self._unsavedLinks.discard(self)

            # Generate a riak object from our model instance
            obj = self.to_riak()
            obj.store()

            if self.countSaveCalls() == 1:
                unsaved = self._unsavedLinks.copy()
                self._unsavedLinks.clear()

                for model in unsaved:
                    model.save()

    def delete(self):
        """Deletes the model instance from Riak. You may override this, but you will want to call this base method to
        perform the actual delete.

        """
        obj = self.to_riak()
        obj.delete()

    def full_clean(self):
        """Performs a full clean. You may override this, but you will want to make sure that you call `clean_fields()`
        _before_ you call `clean()`.

        """
        self.clean_fields()
        self.clean()

    def clean_fields(self):
        """Calls the validators on all fields in this model. Do not overwrite unless you intend to implement this
        functionality yourself. Instead, you should override `clean()`.

        """
        for field in self._fields.values():
            field.validate(self)

    def clean(self):
        """Performs additional model validation. This is a stub method, designed to be overridden by subclasses of
        Model.

        """
        pass

    @classmethod
    def to_python(cls, riakObj):
        """Converts a RiakObject into a Model instance.

        """
        model = cls()

        #FIXME: Uncommenting the following line (which fixes the fact that objects fetched from Riak have no 'pk' set)
        # causes one of the tests to fail, with a reverse link ending up as None instead of an object! ...WUT?
        #model.pk = riakObj.get_key()

        for key, value in riakObj.get_data().iteritems():
            setattr(model, key, value)

        for link in riakObj.get_links():
            print link
            # Add links.
            tag = link.get_tag()
            if tag not in model._links:
                model._links[tag] = ObjectCollection()

            model._links[tag] += link

        for index in riakObj.get_indexes():
            field = index.get_field()
            value = index.get_value()

            if field in model._indexes:
                if not isinstance(model._indexes[field], list):
                    model._indexes[field] = [model._indexes[field]]
                model._indexes[field].append(value)
            else:
                model._indexes[field] = value

        return model

    def to_dict(self):
        """Converts this model instance's content data into a `dict`.

        This does **not** include links or pure secondary indexes.

        """
        values = dict()

        for field in self._fields.values():
            value = field.renderValue(self)

            if value is not None or field.nullable:
                values[field.name] = value

        return values

    def _expandIndexes(self):
        """Get a list of (index, value) tuples for secondary indexes on this object.

        This isn't a dictionary, because the same index may have multiple values on a given object.

        """
        for key, value in self._indexes.iteritems():
            field = self._fields.get(key, None)
            suffix = field.index_suffix if field is not None else "_bin"

            if isinstance(value, (list, tuple)):
                for val in value:
                    yield key + suffix, val
            else:
                yield key + suffix, value

    def _expandLinks(self):
        """Get a list of RiakLink objects for links on this object.

        This isn't a dictionary, because multiple links may have the same tag on a given object.

        """
        for tag, value in self._links.iteritems():
            for link in value.links():
                yield link

    @staticmethod
    def makeRiakObject(bucket, pk, content, indexes=[], links=[]):
        """Creates a RiakObject with the given information.

        bucket: a string; the bucket of the resulting object
        pk: a string; the key of the resulting object
        content: any value that can be converted to JSON (pass None to avoid setting the object's content)
        indexes: an optional list of (index, value) tuples
        links: an optional list of RiakLink objects

        """
        if isinstance(bucket, basestring):
            bucket = rkit.connection.bucket(bucket)

        obj = bucket.new(pk, content)

        # Set indexes
        for index, value in indexes:
            obj.add_index(index, value)

        # Set links
        for link in links:
            obj.add_link(link)

        return obj

    def to_riak(self):
        """Converts this model instance into a RiakObject.

        In order to properly support reverse links, _this function saves the object to Riak_, then continues processing
        the link relationships. The final object returned _does not have its links committed to the database_.

        """
        obj = riak.RiakObject(rkit.connection, self.objects.bucket, self.pk)
        obj.set_data(self.to_dict())

        # Set indexes
        for key, value in self._indexes.iteritems():
            field = self._fields.get(key, None)
            suffix = field.index_suffix if field is not None else "_bin"

            if isinstance(value, (list, tuple)):
                for val in value:
                    obj.add_index(key + suffix, val)
            else:
                obj.add_index(key + suffix, value)

        # Ensure the object is saved, and that we have a primary key
        obj.store()
        self.pk = obj.get_key()

        # Now we can handle links, and set the reverse relationship.

        # Set links
        for key, value in self._links.iteritems():
            #print "to_riak: setting links for field {}: {!r}".format(key, value)
            field = self._fields.get(key, None)
            related_name = field.related_name if field is not None else None
            related_model = field.related_model if field is not None else None

            for link in value.links():
                #print "    link:", link.get_key(), link.get_bucket()
                obj.add_link(link)

                if related_model is None and related_name is None:
                    continue

                bucket = link.get_bucket()
                key = link.get_key()

                # First, try to look it up in our instances registry
                model = Model.instances.get((bucket, key))

                # If related_name is an empty string, we set it to the bucket name, as our default.
                if related_name == "":
                    related_name = self.bucket

                # If not found, get it from the database
                if model is None:
                    link._client = rkit.connection
                    related = link.get()

                    # Set up the relationship in the database.
                    related.add_link(obj, related_name)
                    related.store()

                else:
                    # Set up the reverse relationship.
                    if related_name not in model._links:
                        model._links[related_name] = ObjectCollection()

                    if field.one_to_one:
                        model._links[related_name].clear()

                    # '+=' won't work here because Python interprets it as a __setitem__ when used after a [] expr.
                    model._links[related_name].add(self)

                    self._unsavedLinks.add(model)

        return obj

    def __unicode__(self):
        """Return this model instance as a unicode string.

        """
        if self.pk is not None:
            instance = u'pk={!r}'.format(self.pk)
        else:
            instance = u'no pk'
        return u'<{}.{}({})>'.format(type(self).__module__, type(self).__name__, instance)

    def __str__(self):
        """Return this model instance as an ASCII string.

        """
        if self.pk is not None:
            instance = 'pk={!r}'.format(self.pk)
        else:
            instance = 'no pk'
        return '<{}.{}({})>'.format(type(self).__module__, type(self).__name__, instance)

    def __repr__(self):
        """Return a representation of this model instance as an ASCII string.

        """
        if self.pk is not None:
            instance = 'pk={!r}'.format(self.pk)
        else:
            instance = 'no pk'
        return '<{}.{}({}): {!r}>'.format(type(self).__module__, type(self).__name__, instance, self.to_dict())

#----------------------------------------------------------------------------------------------------------------------
