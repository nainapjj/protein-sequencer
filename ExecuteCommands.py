import subprocess
import os
import tempfile

BASH_FILE = "runBash"
INPUT_FILE = "inputFile"


def __constructFilename(base):
    return tempfile.mkstemp(prefix=base, text=True, dir=os.getcwd())


def runCommandWithInput(cmd, programInput):
    # Construct an input file on the fly
    fileHandle, inputFile = __constructFilename(INPUT_FILE)

    with os.fdopen(fileHandle, 'w') as f:
        f.write(programInput)

    # Construct a bash file on the fly
    fileHandle, bashFile = __constructFilename(BASH_FILE)
    with os.fdopen(fileHandle, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("%s < %s\n" % (cmd, inputFile))

    # Change to executable
    os.chmod(bashFile, 0777)

    # Redirect stderr to STDOUT file in order to get synchronized output
    directory = os.path.join(os.getcwd(), bashFile)
    proc = subprocess.Popen([directory], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, err = proc.communicate()

    # Delete the temp files
    os.remove(bashFile)
    os.remove(inputFile)

    return output
