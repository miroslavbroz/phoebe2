#!/bin/sh

mkdir tmp

for F in * ; do awk '($1>=6500-50-0.01) && ($1<6600+50+0.01)' < $F.vis.dat > tmp/$F ; done

