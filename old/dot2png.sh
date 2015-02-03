#!/bin/bash
set -e

if [ ! `command -v dot` ]
then
    if [ -f "/etc/debian_version" ]
    then
        echo "You do not have graphviz installed, try: sudo apt-get install graphviz"
        exit 1
    else
        echo "You do not have graphviz installed, www.graphviz.org"
        exit 1
    fi
fi


if [ -e "$1.dot" ]; then
    echo "Creating $1.png ..."
    dot -Tpng "$1.dot" -o "$1.png"
    echo "Done"
else
    echo -e "usage: ./dot2png <prefix>\n  and <prefix>.dot will be converted to <prefix>.png"
fi

