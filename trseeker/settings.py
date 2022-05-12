#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Trseeker settings loader.
"""
import os
import yaml
import platform 

SETTINGS_FILENAME = "settings.yaml"
NGRAM_LENGTH = 23
NGRAM_N = 100000000


def load_settings():
    """ Load settings from yaml file.
    @return settings
    """
    file_name = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file_name)[0],
                                 SETTINGS_FILENAME)
    with open(settings_file) as fh:
        settings = yaml.load(fh, Loader=yaml.FullLoader)
    myos = platform.system()
    if myos == "Windows":
        settings["trseeker"]["os"] = "WIN"
    elif myos == "Darwin":
        settings["trseeker"]["os"] = "OSX"
    else: 
        settings["trseeker"]["os"] = "NIX"
    return settings


def save_settings(settings):
    """ Save settings to yaml file.
    @param settings: trseeker settings
    """
    file_name = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file_name)[0],
                                 SETTINGS_FILENAME)
    with open(settings_file, "w") as fh:
        yaml.dump(fh, settings)
