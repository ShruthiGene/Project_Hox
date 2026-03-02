#!/usr/bin/env python3
"""
Fix BED files so start/end (columns 2 and 3) are plain integers.
deepTools computeMatrix fails on scientific notation (e.g. 1.3e+07).
Usage: python fix_bed_coords.py <file.bed> [file2.bed ...]
       python fix_bed_coords.py <dir>   (fix all .bed under dir recursively)
"""
import sys
import os

def fix_line(line):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 3:
        try:
            f[1] = str(int(float(f[1])))
            f[2] = str(int(float(f[2])))
        except ValueError:
            pass  # keep as-is if not numeric
    return "\t".join(f) + "\n"

def fix_file(path):
    with open(path, "r") as f:
        lines = f.readlines()
    new_lines = [fix_line(ln) for ln in lines]
    with open(path, "w") as f:
        f.writelines(new_lines)
    print("Fixed:", path)

def main():
    if not sys.argv[1:]:
        print(__doc__)
        sys.exit(1)
    arg = sys.argv[1]
    if os.path.isfile(arg):
        fix_file(arg)
        for p in sys.argv[2:]:
            if os.path.isfile(p):
                fix_file(p)
        return
    if os.path.isdir(arg):
        for root, _dirs, files in os.walk(arg):
            for name in files:
                if name.endswith(".bed"):
                    fix_file(os.path.join(root, name))
        return
    print("Not a file or directory:", arg)
    sys.exit(1)

if __name__ == "__main__":
    main()
