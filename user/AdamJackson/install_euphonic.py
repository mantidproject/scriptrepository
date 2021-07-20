import importlib
import scipy
import site
import subprocess
import sys
import tempfile

# Update pip, the version that ships with Ubuntu will insist on updating numpy
# to a version that may not be binary-compatible with Mantid
print("Updating Pip")
process = subprocess.run([sys.executable, "-m", "pip", "install",
                          "--user",
                          "pip"],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
print(process.stdout.decode('utf-8'))

with tempfile.TemporaryDirectory() as tmpdirname:
    # We need packaging to inspect version numbers, but not
    # after Euphonic is installed, so install to tempdir
    print("Installing Packaging to temporary directory")
    process = subprocess.run([sys.executable, "-m", "pip", "install",
                              "--prefix", tmpdirname,
                              "packaging"],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    print(process.stdout.decode('utf-8'))
    sys.path.insert(1, glob.glob(tmpdirname + '/lib/*/site-packages')    
    importlib.reload(site)
    globals()['packaging'] = importlib.import_module('packaging')
    from packaging import version

    if version.parse(scipy.version.version) < version.parse('1.0.0'):
        print(f"Scipy version is {scipy.version.version}, updating to 1.0.0")
        process = subprocess.run([sys.executable, "-m", "pip", "install",
                                  "--user",
                                  "scipy==1.0.0"],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT)
        print(process.stdout.decode('utf-8'))

    process = subprocess.run([sys.executable, "-m", "pip", "install",
                              "--user",
                              "euphonic"],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT)
    print(process.stdout.decode('utf-8'))

    globals()['euphonic'] = importlib.import_module('euphonic')
    #print("Please restart Mantid in order to use the installed Euphonic library")
