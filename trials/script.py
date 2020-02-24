import subprocess as sub
import numpy as np
cf = 1
fr = ["250"]
cm = 4
dp = ["5D-2"]
with open('choose_file.dat', 'w') as f1:
    f1.write(str(cf))
with open('file_resolution.dat', 'w') as f2:
    f2.write(" ".join(fr))
with open('choose_method.dat', 'w') as f3:
    f3.write(str(cm))
with open('decomposition_parameters.dat', 'w') as f4:
    f4.write(" ".join(dp))
sub.call(["./main.x"])
    
    
