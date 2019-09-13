#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 03:31:53 2019

@author: rodrigo
"""

from string import ascii_letters
import isotopes as data


fraction = {'H1': 2, 'O16':1}

factor = sum(fraction.values())
fraction.update((key, value/factor) for key, value in fraction.items())

print(fraction)

# for i in fraction:
#     fraction[i] /= factor

def massToAtom(fraction):
    nf = lambda w, mm: w / mm
    # total_ratio = 0
    total_ratio = sum(nf(value, getattr(data, key).mass) for key, value in fraction.items())

    # for i in fraction:
    #     total_ratio += nf(fraction[i], getattr(data, i).mass)


    fraction.update((key, nf(value, getattr(data, key).mass) / total_ratio) for key, value in fraction.items())

    # for i in fraction:
    #     fraction[i] = nf(fraction[i], getattr(data, i).mass) / total_ratio

def atomToMass(fraction):
    wf = lambda n, mm: n * mm
    total_ratio = sum(wf(value, getattr(data, key).mass) for key, value in fraction.items())
    # total_ratio = 0
    # for i in fraction:
    #     total_ratio += wf(fraction[i], getattr(data, i).mass)

    fraction.update((key, wf(value, getattr(data, key).mass) / total_ratio) for key, value in fraction.items())
    # for i in fraction:
    #     fraction[i] = wf(fraction[i], getattr(data, i).mass) / total_ratio


atomToMass(fraction)
print(fraction)

massToAtom(fraction)
print(fraction)
