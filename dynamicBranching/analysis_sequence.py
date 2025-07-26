#!/usr/bin/python

import re

sample_text = '[1, 2];[\'\', 2941.08]'
sample_text_2 = '[2];[1148.74, '']'
sample_text_3 = '[];[2330.53, '']'

def split_line(txt):
  re1='(\\[.*?\\])'	# Square Braces 1
  re2='.*?'	# Non-greedy match on filler
  re3='(\\[.*?\\])'	# Square Braces 2
  rg = re.compile(re1+re2+re3,re.IGNORECASE|re.DOTALL)
  m = rg.search(txt)
  if m:
    sbraces1=m.group(1)
    sbraces2=m.group(2)
    return sbraces1, sbraces2
  return None

def remove_square_brackets(txt):
  try:
    f = txt.strip('[]')
  except:
    f = None
  return f

def split_comma(txt):
  try:
    f = [x.strip() for x in txt.split(',')]
  except:
    f = None
  return f

def string_to_float(txt):
  try:
    f = float(txt)
  except:
    f = None
  return f

def parse_line(txt):
  row1_values = list()
  row2_values = list()
  row1_str, row2_str = split_line(txt)
  x = split_comma(remove_square_brackets(row1_str))
  for y in range(0, len(x)):
    row1_values.append(string_to_float(x[y]))
  x = split_comma(remove_square_brackets(row2_str))
  for y in range(0, len(x)):
    row2_values.append(string_to_float(x[y]))
  return row1_values, row2_values
