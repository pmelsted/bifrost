#!/bin/bash
set -e
if [ -e "$1.dot" ]; then
    echo "Creating $1.png ..."
    dot -Tpng "$1.dot" -o "$1.png"
    echo "Done"
else
    echo -e "usage: ./dot2png <prefix>\n  and <prefix>.dot will be converted to <prefix>.png"
fi

