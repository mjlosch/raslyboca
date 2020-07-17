#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2018-01-01/to/2018-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "165.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "target": "u10.nc",
    'format' : "netcdf",
})
