import os
import subprocess

repeats = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.repeats"
distribute = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.distribute"
executable = "/home/benoit/github/RepeatsCounter/build/RepeatsCounter" 

command = [executable, repeats, distribute]
print("")
print("generateExample will execute:")
print(" ".join(command))
print("")
subprocess.check_call(command)
