#! /bin/bash
#
#cp kmeans_utils.h /usr/bin/gcc/include
#
gcc -c -Wall -I ./libs ./libs/kmeans_utils.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv kmeans_utils.o /lib/kmeans_utils.o
#
echo "Normal end of execution."
