#!/usr/bin/python

"""
set atom type for water boxes
"""


#r=open('30box.xyz','r')
#w=open('30box.xyzt','w')

r=open('xopt.xyz','r')
w=open('xopt.xyzt','w')

for line in range(0,2):
 newline=r.readline()
 w.write(newline)


# TIP3P names
for line in r.readlines():
 if "o" in line:
  newline=line.rstrip()+' OW\n'
  w.write(newline)
 if "h" in line:
  newline=line.rstrip()+' HW\n'
  w.write(newline)
 if "O" in line:
  newline=line.rstrip()+' OW\n'
  w.write(newline)
 if "H" in line:
  newline=line.rstrip()+' HW\n'
  w.write(newline)
