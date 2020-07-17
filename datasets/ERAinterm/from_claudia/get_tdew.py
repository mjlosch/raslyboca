#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2015-01-01/to/2017-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "168.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "target": "tdew.nc",
    'format' : "netcdf",
})
