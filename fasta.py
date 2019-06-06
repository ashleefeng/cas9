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
    
    def __iter__(self):
        return self
        
    def isnext(self):
        if self.prev_line == "":
            return False
        else:
            return True
    
    def next(self):
                
        if self.prev_line == "":
            raise StopIteration
        
        if self.last_ident is None:

            line = self.file.readline()
            # Verify it is a header line
            assert line.startswith (">")

            ident = line.lstrip(">").rstrip('\n')
            # ident = line.split()[1:]
        
        else:
             
            ident = self.last_ident

        sequences = []

        while True:
            line = self.file.readline().rstrip("\r\n")
            if line.startswith(">"):
                self.last_ident = line.lstrip('>')
                break
            elif line == "":
                break
            else:
                sequences.append(line)
        
        self.prev_line = line
        return ident, "".join(sequences)