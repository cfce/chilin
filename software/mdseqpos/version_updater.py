#!/usr/bin/env python
"""
SYNOPSIS: this script is intended to be a hg commit hook-script, i.e. it is 
suppose to be run everytime you commit to the repository.
On every commit, it will write the current changeset information to lib/pkg.cfg
This changeset information is then available for the programs to read.

NOTE: you may need to tweak _conf_file so that it finds the right cfg file
SEE the bottom of this file
"""
# BEGIN: YOU PROBABLY SHOULD LEAVE THIS ALONE
import os
import ConfigParser
import subprocess

try:
    _ver = os.environ['HG_NODE']
except:
    _ver = None

def getHGTip():
    """returns the hg tip changeset"""
    ps = subprocess.Popen(["hg", "tip"], stdout=subprocess.PIPE)
    (out, err) = ps.communicate()
    #the version is the second part of the first line
    first_line = out.split("\n")[0]; 
    return first_line.split(" ")[-1];

def update_config(configFile, ver):
    """Tries to read in the configFile, and updates the current version to ver
    IF configFile dne, then a new configFile is written w/ the version
    """    
    configParser = ConfigParser.ConfigParser()
    
    if len(configParser.read(configFile)) == 0:
        configParser.add_section("pkg_info")

    configParser.set("pkg_info", "hg_changeset", ver)
    f = open(configFile, "w")
    configParser.write(f)
    f.close()
#END: YOU PROBABLY SHOULD LEAVE THIS ALONE

#note this might have to be tweaked for every deployment!
_conf_file = os.path.join(os.getcwd(), "lib", "pkg.cfg")

if __name__ == '__main__':
    if not _ver:
        (ignored,_ver) = getHGTip().split(":")

    update_config(_conf_file, _ver)
