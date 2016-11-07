#!/usr/bin/env bash

VERSION=$(grep '^VERSION' setup.py | sed -E "s/.*'([^']*)'/\1/")
RELEASE_DATE=$(date "+%Y-%m-%d")

echo "[${VERSION}] - ${RELEASE_DATE}"
echo "------------------------------"
echo ""

PREVTAG=$(git describe --tags --abbrev=0 HEAD^)
if [ $? -eq 0 ]
then
    REVRANGE="${PREVTAG}..HEAD"
else
    REVRANGE=""
fi
echo "${REVRANGE}"

# render in reStructured
git log ${REVRANGE} --format='* `%h <https://github.com/tomazc/iCount/commit/%H>`_ %s, %b (%aN)' |
sed -E 's/#([0-9]+)*/\`[#\1\] <https:\/\/github.com\/tomazc\/iCount\/pull\/\1>`_/g'

