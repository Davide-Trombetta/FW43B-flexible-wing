#!/usr/bin/env python3
"""
animationSeries.py - ParaView .series files for the wake animation.

OpenFOAM writes one folder per time step, postProcessing/<function>/<time>/<surface>.vtk,
and ParaView does not group files with the same name in different folders into one
animation. This script writes, for every surface, postProcessing/<function>/<surface>.vtk.series
(the ParaView format for a time series: file names and times); open it in ParaView and press
Play, or File > Save Animation for the frames / a video.

Usage (in the case folder):  python3 tools/animationSeries.py [wakePlanes wakeVortices]
"""
import json, os, sys

def main():
    funcs = sys.argv[1:] or ["wakePlanes", "wakeVortices"]
    for f in funcs:
        root = os.path.join("postProcessing", f)
        if not os.path.isdir(root):
            print(f"{root}: not found, skipped"); continue
        times = []
        for d in os.listdir(root):
            try: times.append((float(d), d))
            except ValueError: pass
        times.sort()
        names = sorted({n for _, d in times for n in os.listdir(os.path.join(root, d)) if n.endswith(".vtk")})
        for n in names:
            files = [{"name": os.path.join(d, n), "time": t} for t, d in times
                     if os.path.isfile(os.path.join(root, d, n))]
            out = os.path.join(root, n + ".series")
            with open(out, "w") as fh:
                json.dump({"file-series-version": "1.0", "files": files}, fh, indent=1)
            print(f"{out}: {len(files)} frames, t = {files[0]['time']:g} ... {files[-1]['time']:g} s")

if __name__ == "__main__":
    main()
