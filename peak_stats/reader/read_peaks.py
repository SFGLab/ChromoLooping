#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: zparteka
"""


class Reader:
    head = None
    lines = None
    path = None

    def __init__(self, filename):
        self.filename = filename
        data = self.read_ascii()
        self.head = data[0]
        self.lines = data[1]
        self.path = filename

    def read_ascii(self):
        f = open(self.filename, "r")
        lines = f.readlines()
        f.close()
        head = lines[0].split("\t")
        lines = lines[1:]
        for i in range(len(lines)):
            lines[i] = lines[i].split()
        return head, lines
