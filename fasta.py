#!/usr/bin/env python

"""
Parse every FASTA record from a file and print.
"""

import sys

class FASTAReader(object):
    
    def __init__(self, file):
        self.file = file
        self.last_ident = None
        self.prev_line = True
        self.isnext = True
    
    def __iter__(self):
        return self
    
    def next(self):
                
        if self.prev_line == "":
            raise StopIteration
        
        if self.last_ident is None:

            line = self.file.readline()
            # Verify it is a header line
            assert line.startswith (">")

            ident = line.lstrip(">").rstrip("\r\n")
            # ident = line.split()[1:]
        
        else:
             
            ident = self.last_ident

        sequences = []

        while True:
            line = self.file.readline()
            if line == '':
                self.isnext = False
                break

            line = line.rstrip("\r\n")
            
            if line.startswith(">"):
                self.last_ident = line.lstrip('>')
                break

            else:
                sequences.append(line)
        
        self.prev_line = line
        return ident, "".join(sequences)