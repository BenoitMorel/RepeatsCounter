import os
import subprocess

msa = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.phy"
part = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.part"
newick = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.newick"
output = "/home/benoit/github/RepeatsCounter/example/small_dataset/small.repeats"
executable = "/home/benoit/github/RepeatsCounter/MSAConverter/build/convert" 

command = [executable, msa, part, newick, output]
print("")
print("generateExample will execute:")
print(" ".join(command))
print("")
subprocess.check_call(command)
