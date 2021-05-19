""" generate_ARGUS_inst_def
Author: James Lord
Generates IDF for ARGUS
"""
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

# generate ARGUS tables
print ('''  <!-- detector components --> 
 <component type="all-rings"  idlist="all">
   <location />
  </component>
  
  <type name="all-rings">
     <component type="rings-front" >
       <location />
     </component>
     <component type="rings-back" >
       <location />     
     </component>    
  </type>''')

for (bank,zsign) in (("front",-1),("back",1)):
    
    print ('  <type name="rings-{bank:s}" >'.format(bank=bank))
    for ring in range(1,7):
        zpos=[None,0.080,0.075,0.070,0.065,0.060,0.055][ring]*zsign
        print ('     <component type="{bank:s}-{ring:d}">'.format(bank=bank,ring=ring))
        print ('<location z="{z:f}" name="ring-{bank:s}-{ring:d}" />'.format(bank=bank,ring=ring,z=zpos))
        print ('</component>')
    print ('</type>')
  

detorder=[]
for bank in ("front","back"):
    for ring in range(1,7):
        r=[None,0.050,0.055,0.060,0.065,0.070,0.075][ring]
        print (' <type name="{bank:s}-{ring:d}">'.format(bank=bank,ring=ring))
        print ('   <component type="ring{ring:d}-pixel">'.format(ring=ring))
        for itheta in range(16):
            if bank=="front":
                theta=((20-itheta)%16)*360.0/16.0
                base=0
            else:
                theta=((itheta+4)%16)*360.0/16.0
                base=96
            det=itheta*6+ring+base
            print ('<location r="{r:f}" t="90" p="{theta:f}" rot="{theta:f}" name="det{det:d}"/>'.format(r=r,theta=theta,det=det))
            detorder.append(det)
        print ('    </component>')
        print ('  </type>')  

for ring in range(1,7):
    halfwidth=[None,0.009,0.010,0.011,0.012,0.013,0.014][ring]
    print(
'''  <type name="ring{ring:d}-pixel" is="detector">  
    <!-- It is implicitely assumed here that the front y-z plane (looking down
         the x-axis) is the surface that see the neutron beam.
         This surface is {width:f}mm along y and 30mm along z and the dept along x is 5mm.  -->
    <cuboid id="ring{ring:d}-pixel-shape">
      <left-front-bottom-point x="0.0" y="-{halfwidth:f}" z="-0.015"  />
      <left-front-top-point  x="0.0" y="{halfwidth:f}" z="-0.015"  />
      <left-back-bottom-point  x="0.005" y="-{halfwidth:f}" z="-0.015"  />
      <right-front-bottom-point  x="0.0" y="-{halfwidth:f}" z="0.015"  />
    </cuboid>
    <algebra val="ring{ring:d}-pixel-shape" />     
  </type>
'''.format(ring=ring,halfwidth=halfwidth,width=2*halfwidth))

print('  <!-- DETECTOR ID LISTS -->')
print (' ')
print ('  <idlist idname="all">')

while detorder:
    n=2
    delta=detorder[1]-detorder[0]
    while n<len(detorder) and detorder[n]==detorder[0]+n*delta:
        n=n+1
    print ('    <id start="{:d}" step="{:d}" end="{:d}"/>'.format(detorder[0],delta,detorder[n-1]))
    detorder=detorder[n:]

print ('  </idlist>')

print ('</instrument>')
