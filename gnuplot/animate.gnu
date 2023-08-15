
# Set terminal to x11

set terminal x11

# Find number of data sets using a perl call, and save it in N.

# IMPORTANT:  YOU NEED TO CHANGE THE NAME OF THE FILE HERE TO
# THE ONE YOU REALLY WANT TO PLOT.

N = `perl -ne '$count += tr/#//; END{print "$count\n"}' alpha.rl`

# Set range in y.

set xrange [0:50]
set yrange [0:1.2]

# Initialize loop variable i.

i=0

# Load looper code.

load "loop.gnu"
