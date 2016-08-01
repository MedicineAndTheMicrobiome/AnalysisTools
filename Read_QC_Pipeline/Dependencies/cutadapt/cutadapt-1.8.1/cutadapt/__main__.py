#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function, division, absolute_import
import sys

try:
	import _preamble
except ImportError:
	pass

from cutadapt.scripts import cutadapt
cutadapt.main()
