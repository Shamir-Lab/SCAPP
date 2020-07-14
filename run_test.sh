#!/usr/bin/env bash

echo "Running test... (May take a few minutes...)"

scapp -g test/test.fastg -b test/test.bam -o test/test_out > /dev/null 2>&1

cd test

python testdiff.py
if [ $? -ne 0 ]
then
  exit 1
fi

rm -r test_out
