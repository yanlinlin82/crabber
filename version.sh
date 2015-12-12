#!/bin/bash

if [ -z "$1" ]; then
	echo "Usage: $0 <version.h.in> [version.h]"
	exit 1
fi

INPUT=$1
OUTPUT=/dev/stdout
if [ -n "$2" ]; then
	OUTPUT=$2
fi

if (git rev-parse 2>/dev/null); then
	# combine number and hash value
	VER="$(git rev-list --all --count)-$(git rev-parse --short HEAD)"
else
	VER="unknown"
fi

SKIP="no"
if [ -f ${OUTPUT} ]; then
	if (diff -q <(cat ${INPUT} | sed 's/\${VER}/'${VER}'/g') ${OUTPUT} >/dev/null); then
		SKIP="yes"
	fi
fi

if [ ${SKIP} == "no" ]; then
	echo "Generating '${OUTPUT}'..." >/dev/stderr
	cat ${INPUT} | sed 's/\${VER}/'${VER}'/g' > ${OUTPUT}
fi

exit 0
