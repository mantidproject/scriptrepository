# === This works inside the MantiPlot shell on Win 7 ===
import subprocess
# install phase
print(subprocess.Popen("python -m pip install --upgrade pip",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate())
print(subprocess.Popen("python -m pip install --user --upgrade mulfc muesr",shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE).communicate())
# use phase
import sys, site
sys.path.insert(0,site.USER_SITE)


# === This works on Mantid Notebook ===
import pip
pip.main("install --user --upgrade muLFC muesr".split())