# -*- coding: utf-8 -*-
"""R-Kit model manager

Copyright (c) 2012-13 Christopher S. Case and David H. Bronke
Licensed under the MIT license; see the LICENSE file for details.

"""
import json

import rkit

from riak import RiakObject


class ModelManager(object):
    def __init__(self, cls):
        self.cls = cls
        self._bucket = None

    @property
    def bucket(self):
        if self._bucket is None:
            self._bucket = rkit.connection.bucket(self.cls.Meta.bucket)

            # Enable searching on this bucket, since we depend on it.
            if not self._bucket.search_enabled():
                self._bucket.enable_search()

        return self._bucket

    def all(self):
        return self.listToPython(self.bucket.get(key) for key in self.bucket.get_keys())

    def get(self, **kwargs):
        result = self.filter(**kwargs)

        if len(result) > 1:
            raise RuntimeError("More than one result from 'get' call.")

        return result[0]

    def filter(self, **kwargs):
        if "pk" in kwargs:
            if len(kwargs) > 1:
                print "Warning! Extra arguments specified with 'pk'. Ignoring."
            return [self.asPython(self.bucket.get(kwargs['pk']))]

        if len(kwargs) == 1:
            field = kwargs.keys()[0]
            value = kwargs.values()[0]

            field = self.cls._fields.get(field)

            if field is None:
                return []

            elif field.searchable:
                return self.listToPython(link.get() for link in self.bucket.search(field.name + ':' + value).run())

            elif field.index:
                return self.listToPython(link.get() for link
                        in rkit.connection.index(self.bucket.get_name(), field.name + field.index_suffix, value).run())

        secondaryIndexes = dict()
        searchable = dict()
        mapreduce = dict()

        # We're filtering more than one thing; time for us to use Map/Reduce or Search
        #query = rkit.connection.add(self.bucket.get_name())
        for f, val in kwargs.iteritems():
            if f == '_search':
                searchable[None] = val

            field = self.cls._fields.get(f)

            if field is None:
                mapreduce[f] = val

            if field.index:
                secondaryIndexes[f + field.index_suffix] = val
            elif field.searchable:
                searchable[f] = val
            else:
                mapreduce[f] = val

        #TODO: Add support for links.

        if searchable:
            query = self.bucket.search(' AND '.join(
                    (field.name + ':' + value) if field.name else value
                    for field, value in searchable
                    ))
        elif secondaryIndexes:
            # Remove the first secondary index from the dictionary, and use that to start the MapReduce query.
            field, val = secondaryIndexes.popitem()
            query = rkit.connection.index(self.bucket.get_name(), field, val)
        else:
            # Start a generic MapReduce query against this bucket.
            query = rkit.connection.add(self.bucket.get_name())

        expr = ['obj.' + ' '.join(self.resolveComparison(field, val)) for field, val in mapreduce.iteritems()]
        expr.extend(('v.metadata.index.' + ' '.join(self.resolveComparison(field, val))
                for field, val in secondaryIndexes.iteritems()))

        mapFunc = """function(v)
                {
                    var matchingValues = v.values.filter(function(val)
                            {
                                var obj = JSON.parse(val.data);
                                return %s;
                            });

                    if(matchingValues.length == 1)
                    {
                        return [{
                            bucket: v.bucket,
                            key: v.key,
                            data: matchingValues[0].data,
                            metadata: matchingValues[0].metadata
                        }];
                    }
                    else
                    {
                        return [];
                    } // end if
                }""" % " && ".join(expr)

        print "Running MapReduce query with map function:", mapFunc
        query.map(mapFunc)
        results = query.run()

        return self.mapredToModels(results)

    def getComparison(self, key):
        if key.endswith('__lt'):
            oper = '<'
        elif key.endswith('__le'):
            oper = '>='
        elif key.endswith('__gt'):
            oper = '>'
        elif key.endswith('__ge'):
            oper = '>='
        elif key.endswith('__ne'):
            oper = '!='
        else:
            return key, '=='

        return key[:-4], oper

    def resolveComparison(self, key, value):
        key, oper = self.getComparison(key)

        return key.replace('__', '.'), oper, json.dumps(value)

    def mapredToModels(self, objects):
        models = list()

        objects = objects or []

        for obj in objects:
            riakObj = self.cls.makeRiakObject(
                    obj['bucket'], obj['key'],
                    None,  # Don't set content here; set the encoded content below.
                    obj['metadata'].get('index', {}).iteritems(),
                    obj['metadata'].get('Links', [])
                    )

            riakObj.set_encoded_data(obj['data'])

            # Add metadata other than index and links. (unneeded?)
            #metadata = riakObj.get_metadata()
            #metadata.update(
            #        (k, v)
            #        for k, v in obj['metadata']
            #        if k not in ('index', 'Links')
            #        )
            #riakObj.set_metadata(metadata)

            models.append(self.asPython(riakObj))

        return models

    def listToPython(self, objects):
        models = list()

        for riakObj in objects:
            models.append(self.asPython(riakObj))

        return models

    def asPython(self, riakObj):
        if riakObj is None or riakObj.get_data() is None:
            return None
        return self.cls.to_python(riakObj)
