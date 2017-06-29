#!/bin/bash

touch src_diff.txt
for file in *.c
do
	touch ${file}diff.txt
	
	diff $file ./binary_src_arnaud_orig/$file >> ${file}diff.txt
	sed -i -e "1s/^/$file \n/" ${file}diff.txt
	sed -i -e '$G' ${file}diff.txt
	cat ${file}diff.txt >>src_diff.txt
	rm ${file}diff.txt
done

for file in *.h
do
	touch ${file}diff.txt
	
	diff $file ./binary_src_arnaud_orig/$file >> ${file}diff.txt
	sed -i -e "1s/^/$file \n/" ${file}diff.txt
	sed -i -e '$G' ${file}diff.txt
	cat ${file}diff.txt >>src_diff.txt
	rm ${file}diff.txt
done