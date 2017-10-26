#!/bin/bash

function showHelp {
version=0.0.1
versionDate="2014-07-07"

echo "$0 - tic/toc timer pair
Usage: $0 [id]              Stores initial time (w/optional id marker)
       $0 [-u|--unique]     Creates and returns unique id
       Using optional ID is recommended, as it allows simulatenous usage.
Notes:
       ttoc [id]  (displays delta time since ttic was called)
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
        -u|--unique)
            shift
            hasId=yes
            isGenerated=yes
            id=$(tr -dc "[:alpha:]" < /dev/urandom | head -c 8)
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
    if [[ $isGenerated ]] ; then
        echo $id
    fi
else
    tmpfile=/tmp/.ttic.time
fi

echo $(($(date +%s%N)/1000000)) > $tmpfile
