import os,sys
from genie_python.genie import *
sys.path.append(r'u:\test\PythonControl')
from PolrefRoutines import *
#############

while True:
    factor=1
    cset(s1vg=1.0*factor, s2vg=0.5*factor, s3vg=1.0*factor, s4vg=1.0*factor)
    waitformove()
    cset(twotheta=1.0*factor, stheta=0.5*factor)
    waitformove()
    pol_run(u=1000, d=1000, total=18000*factor,Title="Ni Fragment in yoke 2th=1.0 dq/q=2.78%", SampleName='Ni Fragment 1')

    factor=2
    cset(s1vg=1.0*factor, s2vg=0.5*factor, s3vg=1.0*factor, s4vg=1.0*factor)
    waitformove()
    cset(twotheta=1.0*factor, stheta=0.5*factor)
    waitformove()
    change(title="Ni Fragment in yoke 2th=2.0 dq/q=2.78%" )
    pol_run(u=1000, d=1000, total=18000*factor)

    factor=4
    cset(s1vg=1.0*factor, s2vg=0.5*factor, s3vg=1.0*factor, s4vg=1.0*factor)
    waitformove()
    cset(twotheta=1.0*factor, stheta=0.5*factor)
    waitformove()
    pol_run(u=1000, d=1000, total=18000*factor,Title="Ni Fragment in yoke 2th=4.0 dq/q=2.78% A",SampleName='Ni Fragment 1')
    pol_run(u=1000, d=1000, total=18000*factor,Title="Ni Fragment in yoke 2th=4.0 dq/q=2.78% B",SampleName='Ni Fragment 1')
