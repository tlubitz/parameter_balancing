# -*- coding: utf-8 -*-

""" Tablib - formats
"""

import _csv as csv
import _json as json
import _xls as xls
import _yaml as yaml
import _tsv as tsv
import _html as html
import _xlsx as xlsx
import _ods as ods
# csv as last, since xlsx and ods might also be detected as csv
available = (json, xls, yaml, tsv, html, xlsx, ods, csv)
