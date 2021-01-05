import os
from datetime import datetime
from pathlib import Path

import numpy as np
import math
import rainflow
import pandas as pd


path_timehistory = Path.cwd()/"timehistoryfiles/results"
files = [f for f in os.listdir(path_timehistory) if not f.startswith('.')]