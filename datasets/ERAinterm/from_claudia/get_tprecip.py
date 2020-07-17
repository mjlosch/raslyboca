#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2017-01-01/to/2019-08-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "228.128",
    "step": "6/12",
    "stream": "oper",
    "time": "00:00:00/12:00:00",
    "type": "fc",
    "target": "total-pre.nc",
    'format' : "netcdf",
})
