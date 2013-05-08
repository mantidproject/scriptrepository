"""
Convert to runs (on Bragg peaks) to d spacing. Use instrument view to select peaks, then the peaks workspace can be used to calculate UB matrix (print_u_v_vectors.py)
"""
instrument='SEQ'
runs=[38629,38630]
rotation_axis='CCR13VRot'
for run in runs:
	runstr='%s_%d'%(instrument,run)
	out_name='d_%s'%runstr
	Load(Filename=runstr,OutputWorkspace=out_name)
	ConvertUnits(InputWorkspace=out_name,OutputWorkspace=out_name,Target='dSpacing')
	axisstr='%s,0,1,0,1' %rotation_axis
	SetGoniometer(Workspace=out_name,Axis0=axisstr)
	Rebin(InputWorkspace=out_name,OutputWorkspace=out_name,Params='0.5,0.005,10',PreserveEvents='0')
