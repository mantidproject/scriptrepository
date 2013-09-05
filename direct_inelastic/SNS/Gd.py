"""
Compact script to reduce Gd data to MD. Used for generating images for releases
"""
import os

run=11499
#initial energy guess
Eguess=50
#Energy bins:    Emin,Estep,Emax
Erange="-10.0,0.2,45.0"
#Data directory	
datadir="/SNS/SEQ/IPTS-4783/data/"
#Output directory
outdir="/SNS/users/3y9/Gd/"
#Run via histogramming
dohist=True
#Do projections
doproj=False

#GetEi and T0 function: note that in different versions of Mantid
#the Monitor1Spec and monitor2Spec might be "0" and "1"
def GetEiT0(ws_name,EiGuess):
    # Run GetEi algorithm
    [Ei, MonPeak,
     MonIndex, Tzero] = GetEi(InputWorkspace=ws_name, Monitor1Spec="1",
                              Monitor2Spec="2", EnergyEstimate=EiGuess)
    # Extract incident energy and T0
    return [Ei,-Tzero]	
	
filename=datadir+"SEQ_"+str(run)+"_event.nxs"
# Load an event Nexus file
LoadEventNexus(Filename=filename,OutputWorkspace="IWS")
# FilterBadPulses(InputWorkspace="IWS",
# OutputWorkspace = "IWS",LowerCutoff = 50)
LoadNexusMonitors(Filename=filename, OutputWorkspace="MonWS")
# Get Ei and -T0 using the function defined before
[Efixed,T0]=GetEiT0("MonWS", Eguess)
# Change all TOF by -T0
ChangeBinOffset(InputWorkspace="IWS", OutputWorkspace="OWS", Offset=T0)
# normalize by proton charge
NormaliseByCurrent(InputWorkspace="OWS",OutputWorkspace="OWS")
# The algorithm for He3 tube efficiency requires wavelength units
ConvertUnits(InputWorkspace="OWS", OutputWorkspace="OWS", 
             Target="Wavelength", EMode="Direct", EFixed=Efixed)	
# Apply correction due to absorption in He3
He3TubeEfficiency(InputWorkspace="OWS", OutputWorkspace="OWS")
# Switch  to energy transfer
ConvertUnits(InputWorkspace="OWS", OutputWorkspace="OWS",
             Target="DeltaE", EMode="Direct", EFixed=Efixed)
# Apply k_i/k_f factor
CorrectKiKf(InputWorkspace="OWS", OutputWorkspace="OWS")
if dohist:
    # Make sure the bins are correct
    Rebin(InputWorkspace="OWS", OutputWorkspace="OWS",
          Params=Erange, PreserveEvents=False)
    # Convert to differential cross section by dividing by the
    # energy bin width
    ConvertToDistribution(Workspace="OWS")
    
# Load vanadium file
LoadNexus(Filename=outdir+"van.nx5", OutputWorkspace="VAN")
# apply overall mask
MaskDetectors(Workspace="OWS", MaskedWorkspace="VAN")
# normalize by Vanadium, if desired
Divide(LHSWorkspace="OWS", RHSWorkspace="VAN",
       OutputWorkspace="OWS")

# Rename workspace to something meaningful
ws_name = "SEQ_"+str(run)
RenameWorkspace(InputWorkspace="OWS",  OutputWorkspace=ws_name)
# Need to fix the goniometer angle by 49.73 degrees
w = mtd[ws_name]
psi = w.getRun()["CCR13VRot"].getStatistics().mean + 49.73
AddSampleLog(Workspace=ws_name, LogName="CCR13VRot_Fixed",
             LogType="Number Series", LogText=str(psi))
# Set the Goiniometer information
SetGoniometer(Workspace=ws_name, Axis0="CCR13VRot_Fixed,0,1,0,1")
# Set the information for the UB matrix
SetUB(Workspace=ws_name,
      a=3.643, b=3.643, c=5.781, alpha=90, beta=90, gamma=120,
      u='1,1,0', v='0,0,1')
# Create the MDEventWorkspace
md_ws_name=ws_name+'_md'
tag=""
if dohist:
    tag += "h"
else:
    tag += "e"
if doproj:
    tag += "wp"
else:
    tag += "np"

md_ws_name += "_" + tag
if not doproj:
    ConvertToMDEvents(InputWorkspace=ws_name, OutputWorkspace=md_ws_name,
                      QDimensions='Q3D', MinValues='-5,-5,-5,-10',
                      QConversionScales='HKL',
                      MaxValues='5,5,5,45', MaxRecursionDepth='1')
else:
    ConvertToMDEvents(InputWorkspace=ws_name, OutputWorkspace=md_ws_name,
                      QDimensions='Q3D', MinValues='-5,-5,-5,-10',
                      QConversionScales='HKL',
                      MaxValues='5,5,5,45', MaxRecursionDepth='1',
                      Uproj='1,1,0', Vproj='1,-1,0', Wproj='0,0,1')
    
# Remove SPE workspace
DeleteWorkspace(Workspace=ws_name)
# Save MDEventWorkspace to a file
md_path = os.path.join(outdir, md_ws_name+'.nxs')
SaveMD(InputWorkspace=md_ws_name, Filename=md_path)
# Delete MDEventWorkspace for safety
#DeleteWorkspace(Workspace=md_ws_name)
logger.notice("Finished run " + str(run))

