"""
Routines used to move packs
"""
#from copy import deepcopy
import numpy as np
from mantid.simpleapi import*
import matplotlib as mpl
import mpl_toolkits.mplot3d as m3d

def V3D2numpy(V3Dobj):
   return np.array([V3Dobj.X(),V3Dobj.Y(),V3Dobj.Z()])

def rot_x(theta):
    M=np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
    return M

def rot_y(theta):
    M=np.array([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
    return M
       
def rotate_pack(wkspc,Dids,cp,ang,relative=False):
   """
   rotate pack at a fixed distance around a vertical axis at the sample position   
   """
   cp_noy= cp-np.array([0,cp[1],0])  #center of pack brought into plane
   ang_from_z=np.arccos(np.dot(cp_noy,np.array([0,0,1]))/np.linalg.norm(cp_noy)) #angle away from z
   rot_m=rot_y(-ang_from_z)
   rot_n=rot_y(ang)
   for idx in Dids:
      h_det=wkspc.getInstrument().getDetector(idx)
      pix_pos=V3D2numpy(h_det.getPos())
      tmppixpos=pix_pos-cp
      if relative==False:
        tmppixpos=np.dot(rot_m,pix_pos)
	cprot=np.dot(rot_m,cp)
      else:
        cprot=cp		
      tmppixpos=np.dot(rot_n,tmppixpos)
      cprot=np.dot(rot_n,cprot)
      tmppixpos=tmppixpos+cprot
      MoveInstrumentComponent(wkspc,DetectorID=idx,X=tmppixpos[0],Y=tmppixpos[1],Z=tmppixpos[2],RelativePosition=False)


def tilt_pack(wkspc,Dids,cp,angle):
   """
   tilt pack about an axis that is along its center in the plane
   """
   cp_noy= cp-np.array([0,cp[1],0])  #center of pack brought into plane
   ang_from_z=np.arccos(np.dot(cp_noy,np.array([0,0,1]))/np.linalg.norm(cp_noy)) #angle away from z
   rot_m=rot_y(-ang_from_z)
   rot_m_b=rot_y(ang_from_z)
   rot_m2=rot_x(-np.radians(angle))
   cp_noy_no_rot=np.dot(rot_m,cp_noy)

   for idx in Dids:
      h_det=wkspc.getInstrument().getDetector(idx)
      pix_pos=V3D2numpy(h_det.getPos())
      tmppixpos=np.dot(rot_m,pix_pos)
      tmppixposa0=tmppixpos-np.dot(rot_m,cp)
      tmppixpos=np.dot(rot_m2,tmppixposa0)
      tmppixpos=tmppixpos+np.dot(rot_m,cp)
      tmppixpos=np.dot(rot_m_b,tmppixpos)
      MoveInstrumentComponent(wkspc,DetectorID=idx,X=tmppixpos[0],Y=tmppixpos[1],Z=tmppixpos[2],RelativePosition=False)
   
    
def calc_pack_center(wkspc,Dids):
   """
   """
   pixpos=np.zeros((3,len(Dids)))
   for idx,Did in enumerate(Dids):
      h_det=wkspc.getInstrument().getDetector(Did)
      pixpos[:,idx]=V3D2numpy(h_det.getPos())
   return pixpos.mean(axis=1)   

def plot_pack(wkspc,Dids):
   pixpos=np.zeros((3,len(Dids)))
   for idx,Did in enumerate(Dids):
      h_det=wkspc.getInstrument().getDetector(Did)
      pixpos[:,idx]=V3D2numpy(h_det.getPos())
   h_f1=mpl.pyplot.figure()
   ax1=m3d.Axes3D(h_f1)   
   ax1.scatter(pixpos[0,:],pixpos[1,:],zs=pixpos[2,:],c='b',marker='o')
   ax1.set_ylabel('y')
   ax1.set_xlabel('x')
   ax1.set_zlabel('z')
   show()
   return pixpos
   
def conv2_hist_dspace(iws,ows,dbin="2.0,0.01,10.0"):
    d_ws2=ConvertUnits(InputWorkspace=iws,Target='dSpacing',EMode='Elastic')
    Rebin(InputWorkspace=d_ws2,OutputWorkspace=ows,Params=dbin,PreserveEvents=False)

def errobar_wkspc(wkspc,specnum,fmtin='bo'):
     """
     plot a single spectrum in matplot lib.
     """    
     x=wkspc.readX(specnum)
     x=(x[1:]+x[:-1])/2.
     y=wkspc.readY(specnum)
     err=wkspc.readE(specnum)
     mpl.pyplot.figure()
     mpl.pyplot.errorbar(x,y,yerr=err,fmt=fmtin)
       
def plot_wkspc(wkspc,specnum,fmtin='bo'):         
     """
     plot a single spectrum in matplot lib.
     """    
     x=wkspc.readX(specnum)
     x=(x[1:]+x[:-1])/2.
     y=wkspc.readY(specnum)
     err=wkspc.readE(specnum)
     mpl.pyplot.figure()
     mpl.pyplot.plot(x,y,fmtin)
     
def spec2numpy(wkspc,specnum):
     x=wkspc.readX(specnum)
     x=(x[1:]+x[:-1])/2.
     y=wkspc.readY(specnum)
     err=wkspc.readE(specnum)
     return (deepcopy(x),deepcopy(y),deepcopy(err))
         
def d_plot(wkspc,s_wid,e_wid):
    d_sum=SumSpectra(InputWorkspace=wkspc,StartWorkspaceIndex=s_wid,EndWorkspaceIndex=e_wid)
    plot_wkspc(d_sum,0)
    
def d_2numpy(wkspc,s_wid,e_wid):
    d_sum=SumSpectra(InputWorkspace=wkspc,StartWorkspaceIndex=s_wid,EndWorkspaceIndex=e_wid)
    out=spec2numpy(d_sum,0)
    return out  
     
