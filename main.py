from script import run
import sys

if len(sys.argv) == 2:
    run(sys.argv[1])
elif len(sys.argv) == 1:
    run(raw_input("enter mdl file name: \n"))
else:
    print "Too many arguments."
