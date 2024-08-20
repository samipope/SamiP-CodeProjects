# -*- coding: utf-8 -*-
"""R-Kit main module

Copyright (c) 2012-13 Christopher S. Case and David H. Bronke
Licensed under the MIT license; see the LICENSE file for details.

"""
import riak
import models

connection = None

#-----------------------------------------------------------------------------------------------------------------------

def connect(host, port=None, transport="PB"):
    """Connects to the riak database specified by `hostname` (_string_) and `port` (_integer_). Optionally, takes
    `transport`. (Valid values are "PB" or "HTTP". Defaults to "PB".)

    """
    global connection

    if port is None:
        if transport == "PB":
            port = 8087
        else:
            port = 8099

    if transport == "PB":
        transport = riak.RiakPbcTransport

    elif transport == "HTTP":
        transport = riak.RiakHttpTransport

    connection = riak.RiakClient(host=host, port=port, transport_class=transport)

#-----------------------------------------------------------------------------------------------------------------------
