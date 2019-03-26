#!/usr/bin/env python
import re
import datetime
import os, hashlib, time, base64, string
from collections import OrderedDict



import pymongo

def connect_to_mongodb(db_obj):

    try:
        client = pymongo.MongoClient('mongodb://localhost:27017',
            username=db_obj["mongodbuser"],
            password=db_obj["mongodbpassword"],
            authSource=db_obj["mongodbname"],
            authMechanism='SCRAM-SHA-1',
            serverSelectionTimeoutMS=10000
        )
        client.server_info()
        dbh = client[db_obj["mongodbname"]]
        return dbh, {}
    except pymongo.errors.ServerSelectionTimeoutError as err:
        return {}, {"taskstatus":0, "errormsg": "open-connection-failed"}
    except pymongo.errors.OperationFailure as err:
        return {}, {"taskstatus":0, "errormsg": "mongodb-auth-failed"}



def order_json_obj(json_obj, ordr_dict):
    
    for k1 in json_obj:
        ordr_dict[k1] = ordr_dict[k1] if k1 in ordr_dict else 1000
        if type(json_obj[k1]) is dict:
            for k2 in json_obj[k1]:
                ordr_dict[k2] = ordr_dict[k2] if k2 in ordr_dict else 1000
                if type(json_obj[k1][k2]) is dict:
                    for k3 in json_obj[k1][k2]:
                        ordr_dict[k3] = ordr_dict[k3] if k3 in ordr_dict else 1000
                    json_obj[k1][k2] = OrderedDict(sorted(json_obj[k1][k2].items(),key=lambda x: float(ordr_dict.get(x[0]))))
                elif type(json_obj[k1][k2]) is list:
                    for j in xrange(0, len(json_obj[k1][k2])):
                        if type(json_obj[k1][k2][j]) is dict:
                            for k3 in json_obj[k1][k2][j]:
                                ordr_dict[k3] = ordr_dict[k3] if k3 in ordr_dict else 1000
                                for kk in json_obj[k1][k2][j].keys():
                                    ordr_dict[kk] = ordr_dict[kk] if kk in ordr_dict else 1000
                                keyList = sorted(json_obj[k1][k2][j].keys(), key=lambda x: float(ordr_dict[x]))
                                json_obj[k1][k2][j] = OrderedDict(sorted(json_obj[k1][k2][j].items(), key=lambda x: float(ordr_dict.get(x[0]))))
            json_obj[k1] = OrderedDict(sorted(json_obj[k1].items(),key=lambda x: float(ordr_dict.get(x[0]))))
    
    return OrderedDict(sorted(json_obj.items(), key=lambda x: float(ordr_dict.get(x[0]))))





