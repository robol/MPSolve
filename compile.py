#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Compile source files in a clean environment. 
#

import tempfile, subprocess, shutil, os, sys

repository = "/home/leonardo/Sviluppo/MPSolve-2.2"

def debug(message, newline = True):
    """Debug information"""
    print "\033[1m=>\033[0m ", message,
    if newline:
        print ""
    sys.stdout.flush()

if __name__ == "__main__":

    # Check if we need some special target
    if (len(sys.argv)) > 1:
        target = sys.argv[1]
    else:
        target = "all"

    # Save informations about working directory
    working_directory = os.getcwd()

    # Create a temorary directory
    temp_folder = tempfile.mkdtemp()
    debug("Created folder %s" % temp_folder)

    # Clone the respository in there
    debug("Entering %s" % temp_folder)
    os.chdir(temp_folder)

    debug("Cloning git repository...", False)
    p = subprocess.Popen(["git", "clone", repository], stdout=subprocess.PIPE,
                         stderr = subprocess.PIPE)
    (output, errors) = p.communicate()
    if p.wait() != 0:
        print "failed"
        debug("git process ended with errors that are reported here:")
        print ""
        print errors
        print ""
    else:
        print "done"
        debug("Compiling...", False)
        os.chdir("MPSolve-2.2")
        p = subprocess.Popen("make", shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (output, errors) = p.communicate()
        if p.wait() != 0:
            print "failed"
            debug("Compilation ended with errors that are reported here:")
            print ""
            print errors
            print ""

        else:
            print "done"
            if target == "check":
                p = subprocess.Popen("make check", shell=True)
                p.wait()
            elif target == "test":
                p = subprocess.Popen(["bash", "./Maketest","all"])
                p.wait()
            
            debug("Dropping you to a shell with the compiled binaries. Press CTRL+D to ")
            debug("release it and delete the compilation directory (with binaries in it)")

            os.system("bash")

    debug("Cleaning %s" % temp_folder)
    shutil.rmtree(temp_folder)
        

    
