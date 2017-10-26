#!/bin/bash

function showHelp {
version=0.0.1
versionDate="2014-07-07"

echo "$0 - tic/toc timer pair
Usage: $0 [id]
       Stores initial time (w/optional id marker)
Notes:
       ttic [id] (used for marking initial time)
       Temporary timer file is stored at /tmp/.ttic.[ID.]time
Example
       # Global timer (not recommended)
       ttic && (do work) && ttoc
       # Using FooBar as id
       ttic FooBar && (do work) && ttoc FooBar
       # Using a randomly generated id
       id=\$(ttic -u) && (do work) && ttoc \$id
Mainted at: https://gist.github.com/swarminglogic/87adb0bd0850d76ba09f
Author:     Roald Fernandez (github@swarminglogic.com)
Version:    $version ($versionDate)
License:    CC-zero (public domain)
"
    exit $1
}


while test $# -gt 0; do
    case "$1" in
        -h|--help)
            showHelp 0
            ;;
        *)
            hasId=yes
            id=$1
            shift
            break
            ;;
    esac
done

if [[ $hasId ]] ; then
    tmpfile=/tmp/.ttic.${id}.time
else
    tmpfile=/tmp/.ttic.time
fi

if [ ! -e $tmpfile ] ; then
    echo "Did not find initalized time file. Run ttic with same id before ttoc!"
    exit 1
fi

tic=`cat $tmpfile`
toc=$(($(date +%s%N)/1000000))
delta=$(($toc - $tic))

LC_NUMERIC=C LC_COLLATE=C
printf '%g\n' $(bc <<< "scale=3; ${delta}/1000")
