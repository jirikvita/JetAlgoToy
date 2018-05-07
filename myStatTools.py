#!/usr/bin/python

from math import sqrt,pow


def addInQuad(x):
    s = 0.
    y = [a*a for a in x]
    s = sqrt(sum(y)) 
    return s
