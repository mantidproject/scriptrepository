import subprocess
# install phase
print(subprocess.Popen("python -m pip install --upgrade pip",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate())
print(subprocess.Popen("python -m pip install --upgrade --user https://github.com/bonfus/muesr/archive/newcif.zip",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate())