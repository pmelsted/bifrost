#!/bin/bash
set -e
if [ -e "$1.dot" ]; then
    dot -Tpng "$1.dot" -o "$1.png"
    echo "Created $1.png"
else
    echo -e "usage: ./dot2png <prefix>\n  and <prefix>.dot will be converted to <prefix>.png"
fi

