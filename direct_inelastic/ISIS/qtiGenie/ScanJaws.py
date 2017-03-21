import math
######################################################################
# STEP 1: In instrument view choose Bragg peak,  then in the Draw menu select an area around it and "Apply and Save" 
# as ROI to Workspace, which will create a workspace called MaskWorkspace, which is called below.

# STEP 2: Run the following script

######################################################################
# user defined values:
jawstart=48    #jaw width for first run
jawstep= -4    # jaw step for each run
run_no=[32715,32716, 32717, 32718, 32719, 32720, 32721, 32722, 32723, 32724, 32725, 32727, 32728]#[range(32715,32726),32727,32728]
jaw_opening = range(jawstart,jawstart+len(run_no)*jawstep,jawstep)
folder='/archive/NDXMERLIN/Instrument/data/cycle_16_5/' # make folder string empty to use Mantid data search path
instname='MER'
dmin=2000  #you can convert to d-spacing (uncomment ConvertUnits line below) and go for a wide range 
dmax=3000  #of d-spacings to be integrated over or you might want to define it more narrowly in TOF units


######################################################################
intensity=[]
for run,width in zip(run_no,jaw_opening):

    fname=instname+str(run)+'.raw'
    print ' processing file ', fname, 
    
    jawWS = Load(Filename=folder+fname,LoadMonitors='Include')
    jawWS = NormaliseByCurrent(InputWorkspace=jawWS)
    jawWS = MaskDetectors(Workspace=jawWS,MaskedWorkspace='MaskWorkspace') # MaskWorkspace comes from selected area save as ROI to MaskWorkspace
    jawWS = SumSpectra(InputWorkspace='jawWS', IncludeMonitors=False)
    #jawWS = ConvertUnits(InputWorkspace=jawWS,  Target='dSpacing')
    jawWS=Integration(InputWorkspace=jawWS,RangeLower=dmin,RangeUpper=dmax)
    print ': Jaw width= ',width,'Integrated Counts=',jawWS.dataY(0)[0]
    intensity.append(jawWS.dataY(0)[0])

plot(jaw_opening,intensity)
xlim(min(jaw_opening),max(jaw_opening))
ylim(0.0,math.ceil(max(intensity)/10.0**(int(math.log10(max(intensity)))))*10**(int(math.log10(max(intensity)))))
xlabel('Jaw Width (mm)')
ylabel('Integrated Counts (arb. units)')
