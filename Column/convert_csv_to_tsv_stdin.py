#!/usr/bin/env python

import csv, sys
csv.writer(sys.stdout, dialect='excel-tab', lineterminator="\n").writerows(csv.reader(sys.stdin))
