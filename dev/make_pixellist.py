#!/usr/bin/env python3

import sys
import itertools

if __name__ == "__main__":
    x1 = int(sys.argv[1])
    x2 = int(sys.argv[2])
    y1 = int(sys.argv[3])
    y2 = int(sys.argv[4])
    out_fn = sys.argv[5]

    lines = ["%d %d" % (x,y) for x,y in itertools.product(range(x1,x2+1), range(y1,y2+1))]

    with open(out_fn, "w") as of:
        of.write("\n".join(lines))
