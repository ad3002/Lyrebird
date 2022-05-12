#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 13.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:
    
- MongoDBReader(object)
        
"""
try:
    import pymongo
    from pymongo import Connection
except:
    print("WARNING: Install pymongo")
try:
    from bson.code import Code
except:
    print("WARNING: Install bson")


class MongoDBReader(object):
	''' Wrapper for mongodb.'''
	
	def __init__(self):
		pass

	def get_trdb_conn(self):
		''' Connect to **trdb** database.'''
		connection = Connection()
		db = connection.trdb
		return db
