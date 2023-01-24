import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request


run(session, "open 1lz1")
run(session, "open 2hed")
run(session, "torsion #1/A:1@n,ca,cb,cg")
print(run(session, "torsion #1/A:1@n,ca,cb,cg"))
print(run(session, "angle #1/A:1@n #1/A:1@ca #1/A:1@cb"))
