"""
read compressed (.gz .bz) files
"""
#!/usr/bin/env python
# encoding: utf-8
import sys
import bz2
import gzip
import urllib
from subprocess import Popen
from subprocess import PIPE


import sys

if sys.version_info[0] == 3:
	from urllib.request import urlopen
else:
	# Not Python 3 - today, it is most likely to be Python 2
	# But note that this might need an update when Python 4
	# might be around one day
	from urllib import urlopen


def nopen(f, mode="rb"):
	if not isinstance(f, str):
		return f
	if f.startswith("|"):
		p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
		if mode[0] == "r": return p.stdout
		return p
	return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
		else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
		else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
		else urlopen(f) if f.startswith(("http://", "https://","ftp://")) \
		else open(f, mode)

 
def reader(fname):
	for l in nopen(fname):
		yield l.decode('utf8').strip().replace("\r", "")
	
