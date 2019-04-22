#!/bin/bash
#Pulls a quarter of the reads from a fastq file randomly.

getpos () {
	pos=$(((1 + $RANDOM % $inputlength) * 4))
}

inputlength=$(wc -l < $1)
inputlength=$(expr $inputlength / 4)

for x in $(seq $inputlength); do
	getpos
	for y in $(seq 4); do
		$(sed "$($pos, $((pos + $y)))!d" $1) >> "short_$1"
	done
done
