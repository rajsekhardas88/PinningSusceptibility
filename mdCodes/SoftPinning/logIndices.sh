#!/bin/bash

nb=$1 #nblocks
a=$2 #base
r=$3 #runlength

awk "BEGIN{
bs=$r/$nb;
for(i=0;i<$nb;i++){
  printf(\"%d\n\",i*bs)
  x=1
  x0=0
  jmax=log($r-i*bs)/log($a)
  for(j=0;j<=jmax;j++){
    if (x*$a-(x*$a)%1==x0) x+=1; else {x*=$a;x=x-x%1}
    y=x+i*bs
    if (y<$r) printf(\"%d\n\",x+i*bs)
    x0=x-x%1
  }
}
}" | sort -n | awk '{if ($1>x||NR==1) {print $0;x=$1}}'
