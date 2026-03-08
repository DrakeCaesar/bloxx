#!/bin/bash
# Reads all baseline files, counts type transitions, prints a 6x6 table.
# Types: 0=Empty, 1=Blue, 2=Red, 3=Green, 4=Yellow, 5=White

cd "$(dirname "$0")"

awk '
# Find first "(single_digit)" in string s — skips "(r,c)" coords (they have a comma).
function findType(s,    p, c) {
    for (p = 1; p <= length(s); p++) {
        if (substr(s, p, 1) == "(") {
            c = substr(s, p+1, 1)
            if (c ~ /[0-9]/ && substr(s, p+2, 1) == ")") {
                return c + 0
            }
        }
    }
    return -1
}

BEGIN {
    split("Empty,Blue,Red,Green,Yellow,White", names, ",")
    split("empty,blue,red,green,yellow,white", configNames, ",")
    for (i = 1; i <= 6; i++) {
        n = configNames[i]
        while (length(n) < 7) n = n " "
        configNames[i] = n
    }
}

{
    line = $0

    # Strip ANSI escape codes
    while (match(line, /\033\[[0-9;]*m/))
        line = substr(line, 1, RSTART-1) substr(line, RSTART+RLENGTH)

    # Only process Step lines containing a transition
    if (line !~ /^Step [0-9]+:/ || line !~ / -> /) next

    # Split at first " -> "
    arrowPos = index(line, " -> ")
    left  = substr(line, 1, arrowPos - 1)
    right = substr(line, arrowPos + 4)

    # From type: "empty" at end of left means type 0, otherwise find "(N)"
    if (left ~ /empty[[:space:]]*$/) {
        from = 0
    } else {
        from = findType(left)
    }

    # To type: first "(N)" in right half
    to = findType(right)

    if (from >= 0 && to >= 0)
        cnt[from, to]++
}

END {
    # Forward format: from\to (same as config.h FORWARD_TRANSITIONS)
    print "Forward (from -> to)"
    print ""
    printf "%-8s", "from\\to"
    for (j = 0; j <= 5; j++) printf "%8s", names[j+1]
    print ""
    printf "%s\n", "------------------------------------------------------------------------"
    for (i = 0; i <= 5; i++) {
        printf "%-8s", names[i+1]
        for (j = 0; j <= 5; j++)
            printf "%8s", (cnt[i,j]+0 > 0 ? cnt[i,j] : ".")
        print ""
    }
    printf "%s\n", "------------------------------------------------------------------------"
    printf "%-8s", "TOTAL"
    for (j = 0; j <= 5; j++) {
        col = 0
        for (i = 0; i <= 5; i++) col += cnt[i,j]+0
        printf "%8s", (col > 0 ? col : ".")
    }
    print ""

    # Backward format: same 6x6 layout as config.h BACKWARD_TRANSITIONS
    # Backward from i to j = forward from j to i, so cell [i][j] = cnt[j][i]
    print ""
    print "Backward (from -> to)"
    print ""
    printf "%-8s", "from\\to"
    for (j = 0; j <= 5; j++) printf "%8s", names[j+1]
    print ""
    printf "%s\n", "------------------------------------------------------------------------"
    for (i = 0; i <= 5; i++) {
        printf "%-8s", names[i+1]
        for (j = 0; j <= 5; j++)
            printf "%8s", (cnt[j,i]+0 > 0 ? cnt[j,i] : ".")
        print ""
    }
    printf "%s\n", "------------------------------------------------------------------------"
    printf "%-8s", "TOTAL"
    for (j = 0; j <= 5; j++) {
        col = 0
        for (i = 0; i <= 5; i++) col += cnt[j,i]+0
        printf "%8s", (col > 0 ? col : ".")
    }
    print ""

    # Boolean matrices for pasting into config.h (lowercase, no space after {)
    print ""
    print "--- Paste into config.h: FORWARD_TRANSITIONS ---"
    print ""
    for (i = 0; i <= 5; i++) {
        row = "    /* " configNames[i+1] "*/ {"
        for (j = 0; j <= 5; j++) {
            row = row (j > 0 ? ", " : "") ((cnt[i,j]+0) > 0 ? 1 : 0)
        }
        print row "},"
    }
    print ""
    print "--- Paste into config.h: BACKWARD_TRANSITIONS ---"
    print ""
    for (i = 0; i <= 5; i++) {
        row = "    /* " configNames[i+1] "*/ {"
        for (j = 0; j <= 5; j++) {
            row = row (j > 0 ? ", " : "") ((cnt[j,i]+0) > 0 ? 1 : 0)
        }
        print row "},"
    }
}
' level*_type*_baseline.txt
