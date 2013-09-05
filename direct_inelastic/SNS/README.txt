The scripts in this folder are in use (in some form) at SNS direct geometry inelastic instruments (ARCS, CNCS, HYSPEC, SEQUOIA)

Sripts to deal with sample logs and/or experiment logging:

Scripts that split data:  
  new_scan.py - Updated version of rocking curve scan.py, using FilterEvents
  powderTdep_OP.py - Powder line intensity - temperature dependence
  reduce_mult_T.py - Reduce temperature scan. Going as far in events as possible, split, then extract useful information
  scan.py - rocking scan using continuous rotations
  Temp_plot_w.py - temperature scan on a peak using a mask
  T_scan.py - Quick temperature scan
  WB_Tdep.pyTemperature dependence of scattering intensity between Qmin,Qmax

Scripts to do some reduction/analysis:
  AddAndReduce.py - Simple script to add runs together, reduce using DGSReduction, and save to NXSPE format.
  Gd.py - Compact script to reduce Gd data to MD. Used for generating images for releases. Heavy use of workspace groups
  reduce_from_autoreduce_script.py - Reduction including generation of an experiment log 
  reducenoguess.py - Simplified version of reduction using DGSReduction algorithm
  sum_split_plot_slice.py - continuous rotation to a 2D slice in SliceViewer
  Tdep50meV.py - Plot temperature dependence for data in a set of runs
  whitebeampeaks.py - gets UB from a white beam run, ignore those from Aluminum

Scripts to do some postprocessing:
  groupNxspe.py - Reduce file size of NXSPE files by grouping pixels togeter. 
  load_n_plot.py - Quicly look at powder data in NXSPE files
  merge.py - Merge all MD nexus files in a folder
  writeMDhisto.py - Saves an md histo workspace with wname (2D slice from SliceViewer) into an ascii file

Other: 
  average_Ei_T0.py - Print incident energy and t0. Useful to check statistics of monitors in individual runs.
  fixDASLogs.py - function to change initial time for logs to coincide with the first proton pulse
  generate_string_base.py - Generate a list of runs where some log values obey certain criteria. Input can be from log.py 
  log.py - Quicly check logs for multiple runs 
  make_d_spacing.py - Convert to runs (on Bragg peaks) to d spacing. Use instrument view to select peaks, then the peaks workspace can be used to calculate UB matrix (print_u_v_vectors.py)
  MaskBTP.py - Mask banks, tubes, and pixels on ARCS/SEQUOIA/CNCS/HYSPEC 
  print u_v_vectors.py - print u and v vectors for mslice/horace. Use in conjunction with make_d_spacing.py
  testtib.py - Test if time independent background range is OK. Uses MantidPlot plotting
  tof_plot_make.py - Make some plots using matplotlib
  Vec_move.py - Routines used to move packs


