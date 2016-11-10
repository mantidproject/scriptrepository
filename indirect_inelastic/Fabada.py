import os
import subprocess

from mantid.api import *
from mantid.kernel import *
 
class Fabada(PythonAlgorithm):
	def PyInit(self):

		self.declareProperty("Generate_Pfile", False,"Wether to generate a Pfin.par")
		self.declareProperty('Nsteps',1000)
		self.declareProperty('Nwrite',100)
		self.declareProperty('Xmin',-0.5)
		self.declareProperty('Xmax',0.5)
		self.declareProperty('SPCmin',1)
		self.declareProperty('SPCmax',100)
		self.declareProperty('W1',0.0)
		self.declareProperty('W2',0.0)
		self.declareProperty('W3',0.0)
		self.declareProperty('W4',0.0)
		self.declareProperty('W5',0.0)
#		self.declareProperty('Xconv',0.0)
		self.declareProperty('Nlorentz',2,"If the model is SumLorentz, the nuber of Lorentzians")
		self.declareProperty("Model","SumLorentz",StringListValidator(["DifSuml","DifRot"]),"Model to fit the data")
#		self.declareProperty("Background", False,"Wether to substract background")
		self.declareProperty("Delta", False,"Wether to use a delta function")
		self.declareProperty(MatrixWorkspaceProperty("WSFout", "FunctionsOut", direction = Direction.Output), "Workspace with output parameters")
		self.declareProperty(MatrixWorkspaceProperty("Data", "SQ_Clada_T260", direction = Direction.Input), "Workspace with S(q,w) from Data")
		self.declareProperty(MatrixWorkspaceProperty("Resolution", "SQ_Vanadium", direction = Direction.Input), "Workspace with S(q,w) from Vanadium")
		
		

		
	def PyExec(self):
		
		self.log().notice('')
		self.log().notice('###### FABADA ######')
		
		directory = os.path.normpath(config['defaultsave.directory'])
		if not os.path.exists(os.path.join(directory,"fabada.exe")):
			raise RuntimeError("Unable to find Fabada executable in '%s'. Please install executable here." % directory)
		
# getting data to generate the parameter file
		Nlorentz = self.getProperty("Nlorentz").value
		Nsteps = self.getProperty("Nsteps").value
		Nwrite = self.getProperty("Nwrite").value
		Xmin = self.getProperty("Xmin").value
		Xmax = self.getProperty("Xmax").value
		W1 = self.getProperty("W1").value
		W2 = self.getProperty("W2").value
		W3 = self.getProperty("W3").value
		W4 = self.getProperty("W4").value
		W5 = self.getProperty("W5").value
		SPCmin = self.getProperty("SPCmin").value
		SPCmax = self.getProperty("SPCmax").value
#		Ymaxv = self.getProperty("Xconv").value
		Ymaxv = -0.01
		Function = self.getProperty("Model").value
#		Background = self.getProperty("Background").value
		Delta = self.getProperty("Delta").value
		Generate = self.getProperty("Generate_Pfile").value
#		self.log().notice('I will use: '+Function)
		Background = False


			
		
		
# getting the data
		sqw = self.getProperty("Data").value
		NH = sqw.getNumberHistograms()
		NB = sqw.getNumberBins()
#		self.log().notice('Data, NH: '+str(NH)+' NB: '+str(NB))
# GETING q values
		a=sqw.getAxis(1)
		qvalues = a.extractValues()		
#		self.log().notice('qvalues : ' + str(qvalues))

# getting the resolution
		sqwR = self.getProperty("Resolution").value
		NHR = sqwR.getNumberHistograms()
		NBR = sqwR.getNumberBins()
#		self.log().notice('Resolution, NH: '+str(NHR)+' NB: '+str(NBR))
		if NH!=NHR:
			self.log().notice('Puit! I can not fit different number of histograms of data and resolution!!, NH: '+str(NH)+' NHR: '+str(NHR))
		if NB!=NBR:
			self.log().notice('Puit! I can not fit different number of bins of data and resolution!!, NB: '+str(NB)+' NBR: '+str(NBR))
		NH = SPCmax-SPCmin+1
		NHR = SPCmax-SPCmin+1
# Defining the output
# Outputfunctions
		WSFout = WorkspaceFactory.create("Workspace2D",NVectors=NH,XLength=NB,YLength=NB)
		self.setProperty("WSFout", WSFout) 

# definining values for  parameter file

		Nfunc = NH
		Xmaxv=abs(min(Xmin,Xmax))

		outfile = os.path.join(directory, "out.dat")
		chifile = os.path.join(directory, "chi.dat")
		nomMMchain = os.path.join(directory, "MMchain.dat")
		temperature = 1.0
		sigmaexp = 0.0
		Nstepswrite = Nwrite
		varpar = 0
		gpar = 3
		incgpar = 10
		incgchi = 0.5
		prior = 0
		annealing = 0
		stepannealing = 0
		Tini = 10.0
		Tfin = 1.0
		incT = 0.01
		acceptance = 0.33
		Nequal = 0
		Nrel = 0
		Ncont = 0
		Npeak = 0
		
# *************** defining values for parameter file ******************
# *************** defining values for parameter file ******************
# *************** defining values for parameter file ******************
# *************** defining values for parameter file ******************
# *************** defining values for parameter file ******************

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -----------------------------Difussion (or delta) +Lorentzians           -----------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

		if Function == 'DifSuml':

			if Nfunc > 20:
				self.log().notice('I can not fit more than 20 functions, sorry!')
			ifunc=111
			Npar=130
			Parname = ["Y0","P1","P2","W1","W2","W3","W4","W5","diff",
			"A1f1","A1f2","A1f3","A1f4","A1f5","A1f6","A1f7","A1f8","A1f9","A1f10","A1f11","A1f12","A1f13","A1f14","A1f15","A1f16","A1f17","A1f18","A1f19","A1f20",
			"A2f1","A2f2","A2f3","A2f4","A2f5","A2f6","A2f7","A2f8","A2f9","A2f10","A2f11","A2f12","A2f13","A2f14","A2f15","A2f16","A2f17","A2f18","A2f19","A2f20",
			"A3f1","A3f2","A3f3","A3f4","A3f5","A3f6","A3f7","A3f8","A3f9","A3f10","A3f11","A3f12","A3f13","A3f14","A3f15","A3f16","A3f17","A3f18","A3f19","A3f20",
			"A4f1","A4f2","A4f3","A4f4","A4f5","A4f6","A4f7","A4f8","A4f9","A4f10","A4f11","A4f12","A4f13","A4f14","A4f15","A4f16","A4f17","A4f18","A4f19","A4f20",
			"A5f1","A5f2","A5f3","A5f4","A5f5","A5f6","A5f7","A5f8","A5f9","A5f10","A5f11","A5f12","A5f13","A5f14","A5f15","A5f16","A5f17","A5f18","A5f19","A5f20",
			"c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19","c20",
			"tau","alfa","delta"]

# ############ I GUESS THE STARTING PARAMETERS ######
			deltax=abs(sqw.dataX(0)[2]-sqw.dataX(0)[1])
			fmax = max(sqw.dataY(0))
			mid = [i for i in sqw.dataY(0) if i>fmax/4 ]
			yesno = 	 sqw.dataY(0)[50] in mid 
			xtest = []
			for i in range(len(sqw.dataY(0))):
				yesno = 	 sqw.dataY(0)[i] in mid 				
				if yesno == True:
					xtest.append(sqw.dataX(0)[i])
			
			Wg = abs(max(xtest)-min(xtest))/10
			Wguess = [W1,W2,W3,W4,W5]
			if W1 == 0:
				Wguess[0]=Wg
			if W2 == 0:
				Wguess[1]=Wg
			if W3 == 0:
				Wguess[2]=Wg
			if W4 == 0:
				Wguess[3]=Wg
			if W5 == 0:
				Wguess[4]=Wg

			if Delta==True:
				Nlorentz +=1

			if Background == True:
				Par = [0.0,0.0,0.0]
				Parjump = [0.1,0.1,0.0]
				Parlog = ["al","al","al"]
				Parmin = [-1000,-1000,-1000]
				Parmax =[1000,1000,1000]
			if Background == False:
				Par = [0.0,0.0,0.0]
				Parjump = [0.0,0.0,0.0]
				Parlog = ["n","n","n"]
				Parmin = [-1000,-1000,-1000]
				Parmax =[1000,1000,1000]
# I go for the W of lorentzians
			for i in range(Nlorentz):
				Par.append(Wguess[i-1])
				Parjump.append(0.1)
				Parlog.append("l")
				Parmin.append(deltax/2.0)
				Parmax.append(abs(max(abs(Xmin),abs(Xmax)))*10)
				
			for i in range(5-Nlorentz):
				Par.append(6.66)
				Parjump.append(0.0)
				Parlog.append("n")
				Parmin.append(0.0)
				Parmax.append(666.0)
				

# I go fo the diffusion
			if Delta == False:
				Dguess = Wguess/qvalues[0]**2
				Par.append(Dguess)
				Parjump.append(Dguess*0.1)
				Parlog.append("al")
				Parmin.append(0.0)
				Parmax.append(100)
				self.log().notice("My first guess for D is: "+str(Dguess))
			if Delta == True:
				Par.append(0.0)
				Parjump.append(0.0)
				Parlog.append("n")
				Parmin.append(0.0)
				Parmax.append(666)
						
# WE WRITE a's different from zero
			for ilorentz in range(Nlorentz):
				for i in range(Nfunc):
					Par.append(fmax)
					Parjump.append(fmax*0.2)
					Parlog.append('al')
					Parmin.append(fmax/100)
					Parmax.append(fmax*100)
				for i in range(20-Nfunc):
					Par.append(0.0)
					Parjump.append(0.0)
					Parlog.append('n')
					Parmin.append(0.0)
					Parmax.append(fmax*100)
# WE WRITE a's different equal to  zero
			for ilorentz in range(5-Nlorentz):
				for i in range(Nfunc):
					Par.append(0)
					Parjump.append(0)
					Parlog.append('n')
					Parmin.append(0)
					Parmax.append(666)
				for i in range(20-Nfunc):
					Par.append(0.0)
					Parjump.append(0.0)
					Parlog.append('n')
					Parmin.append(0.0)
					Parmax.append(fmax*100)


# WE WRITE C
			for i in range(Nfunc):
				Par.append(0)
				Parjump.append(0.0)
				Parlog.append('al')
				Parmin.append(-10.0)
				Parmax.append(10)

			for i in range(20-Nfunc):
				Par.append(0.0)
				Parjump.append(0.0)
				Parlog.append('n')
				Parmin.append(-6.0)
				Parmax.append(6.0)

# we add delta of x
			Par.append(sqw.dataX(0)[2]-sqw.dataX(0)[1])
			Parjump.append(0.0)
			Parlog.append('n')
			Parmin.append(-666.0)
			Parmax.append(666.0)
#			self.log().notice(str(len(Par)))

			if Delta == True:
				Par[3]=0.0
				Parjump[3]=0.0
				Parmin[3]=-0.666
#				for i in range(9,14):
#					Par[i]=0.0
#					Parjump[i]=0.0
#					Parmin[i]=0.0
				
				

# ----------------------------------------------------------END--------------------------------------------------------		

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -----------------------------Difussion (or delta) + Rotation          -----------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

		elif Function=='DifRot':
			if Nfunc > 20:
				self.log().notice('I can not fit more than 20 functions, sorry!')
			ifunc=112
			Npar=49
			Parname = ["Y0","P1","P2","Rgir","Dr","diff",
			"F1","F2","F3","F4","F5","F6","F7","F8","F9","F10",
			"F11","F12","F13","F14","F15","F16","F17","F18","F19","F20",
			"c1","c2","c3","c4","c5","c6","c7","c8","c9","c10",
			"c11","c12","c13","c14","c15","c16","c17","c18","c19","c20",
			"tau","alfa","delta"]

# ############ I GUESS THE STARTING PARAMETERS ######
			fmax = max(sqw.dataY(0))
			mid = [i for i in sqw.dataY(0) if i>fmax/4 ]
			yesno = 	 sqw.dataY(0)[50] in mid 
			xtest = []
			for i in range(len(sqw.dataY(0))):
				yesno = 	 sqw.dataY(0)[i] in mid 				
				if yesno == True:
					xtest.append(sqw.dataX(0)[i])
				
			Wguess = abs(max(xtest)-min(xtest))/10
			if Background == True:
				Par = [0.0,0.0,0.0]
				Parjump = [0.1,0.1,0.0]
				Parlog = ["al","al","al"]
				Parmin = [-1000,-1000,-1000]
				Parmax =[1000,1000,1000]
			if Background == False:
				Par = [0.0,0.0,0.0]
				Parjump = [0.0,0.0,0.0]
				Parlog = ["n","n","n"]
				Parmin = [-1000,-1000,-1000]
				Parmax =[1000,1000,1000]
# I go fo the radius of giration
			Par.append(1.0)
			Parjump.append(0.1)
			Parlog.append("al")
			Parmin.append(0.0)
			Parmax.append(100)
				
# I go fo the rotation difusion
			Par.append(1.0)
			Parjump.append(0.1)
			Parlog.append("al")
			Parmin.append(0.0)
			Parmax.append(100)
# I go fo the diffusion
			if Delta == False:
				Dguess = Wguess/qvalues[0]**2
				Par.append(Dguess)
				Parjump.append(Dguess*0.1)
				Parlog.append("al")
				Parmin.append(0.0)
				Parmax.append(100)
				self.log().notice("My first guess for D is: "+str(Dguess))
			if Delta == True:
				Par.append(0.0)
				Parjump.append(0.0)
				Parlog.append("n")
				Parmin.append(0.0)
				Parmax.append(666)
						
# WE WRITE F
			for i in range(Nfunc):
				Par.append(fmax)
				Parjump.append(fmax*0.2)
				Parlog.append('al')
				Parmin.append(sqw.dataX(0)[2]-sqw.dataX(0)[1])
				Parmax.append(fmax*100)
			for i in range(20-Nfunc):
				Par.append(0.0)
				Parjump.append(0.0)
				Parlog.append('n')
				Parmin.append(0)
				Parmax.append(fmax*100)


# WE WRITE C
			for i in range(Nfunc):
				Par.append(0)
				Parjump.append(.01)
				Parlog.append('al')
				Parmin.append(-10.0)
				Parmax.append(10)

			for i in range(20-Nfunc):
				Par.append(0.0)
				Parjump.append(0.0)
				Parlog.append('n')
				Parmin.append(-6.0)
				Parmax.append(6.0)
# we add tau			
			Par.append(1.0)
			Parjump.append(0.0)
			Parlog.append('n')
			Parmin.append(-6.0)
			Parmax.append(6.0)
# we add alfa
			Par.append(0.0)
			Parjump.append(0.0)
			Parlog.append('n')
			Parmin.append(-6.0)
			Parmax.append(6.0)
# we add delta of x
			Par.append(sqw.dataX(0)[2]-sqw.dataX(0)[1])
			Parjump.append(0.0)
			Parlog.append('n')
			Parmin.append(-666.0)
			Parmax.append(666.0)

		else:
			self.log().notice('Function '+Function+' does not exist, pech gehabt!')
			sys.exit()
				
# ----------------------------------------------------------END--------------------------------------------------------	


# ****************************** WRTITING PARAMETER FILE *************************************
# ****************************** WRTITING PARAMETER FILE *************************************
# ****************************** WRTITING PARAMETER FILE *************************************
# ****************************** WRTITING PARAMETER FILE *************************************
# ****************************** WRTITING PARAMETER FILE *************************************
		
		if Generate == True:
# I generate the par list
			Parlist = []
			Parlist.append("title "+" 0 1  "+str(Nfunc))
			Parlist.append("nomin   opt  Nint nomout    xmin       xmax    Function")	
		
			for ispc in range (NH):
				outname="f"+str(ispc)+".dat"
				outnameV="V"+str(ispc)+".dat"
				outnameR="R"+str(ispc)+".dat"
				fullpath = os.path.join(directory, outname)
				fullpathV = os.path.join(directory, outnameV)
				fullpathR = os.path.join(directory, outnameR)
				Parlist.append(fullpath+"  cf  0  "+fullpathR+"\t"+str(Xmin)+"\t"+str(Xmax)+"\t"+str(ifunc)) 
				Parlist.append(fullpathV+"\t"+str(Xmaxv)+"\t"+str(Ymaxv))
			Parlist.append(outfile + ' ' + chifile+' '+nomMMchain)
			Parlist.append(str(temperature))
			Parlist.append(str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp)+" "+str(sigmaexp))
			Parlist.append(str(Npar))
			Parlist.append(str(Nsteps) + ' ' + str(Nstepswrite))
			Parlist.append(str(varpar) + ' ' + str(gpar) + ' ' + str(incgpar) + ' ' + str(incgchi) + ' ' + str(prior) )
			Parlist.append(str(annealing) )
			Parlist.append(str(stepannealing) + ' ' + str(Tini) + ' ' + str(Tfin) + ' ' + str(incT) )
			Parlist.append(str(acceptance) )
			Parlist.append ("#i   Par        o  Par            Parjump     Parmin      Parmax")

			for ipar in range(Npar):
				Parlist.append(str(ipar+1) + '\t' + Parname[ipar] + '\t' + Parlog[ipar] + '\t' + str(Par[ipar]) + '\t\t' + str(Parjump[ipar]) + '\t\t'+ str(Parmin[ipar]) + '\t\t' + str(Parmax[ipar]))
			Parlist.append(str(Nequal) )
			Parlist.append(str(Nrel) )
			Parlist.append(str(Ncont) )
			Parlist.append(str(Npeak) )
# I save the par list
			
			fullpath = os.path.join(directory, 'Pfin.par')
			Pfile = open(fullpath,'w')
			for i in range(len(Parlist)):
				Pfile.write(str(Parlist[i])+"\n")
			Pfile.close()
			
			
#		if Generate == False:
#	                ParTable = mtd['ParTable']	
 #
#			fullpath = os.path.join(directory, 'Pfin.par')
#			Pfile = open(fullpath,'r')
#			write = True
#			Parlist = []
#			while True:
#				line = Pfile.readline()
#					
#				if write == True:
#					Parlist.append(line)
#					
#				if "#" in line:
#					break			
#				
#				if line == "":
#					break					
#				
#			Pfile.close()
			
#			for row in ParTable:
#				values = row.values()
#				Parlist.append(str(row['Parjump'])+"\n")
#			for i in range(Npar):
#				Parlist.append(str(ParTable.cell(i,0))+"\t"+str(ParTable.cell(i,1))+"\t"+str(ParTable.cell(i,2))+"\t"+str(ParTable.cell(i,3))+"\t"+str(ParTable.cell(i,4))+"\t"+str(ParTable.cell(i,5))+"\t"+str(ParTable.cell(i,6))+"\n")
#			Parlist.append("0 \n")
#			Parlist.append("0 \n")
#			Parlist.append("0 \n")
#			Parlist.append("0 \n")
#			
#			fullpath = os.path.join(directory, 'Pfino.par')
#			Pfile = open(fullpath,'w')
#			for i in range(len(Parlist)):
#				Pfile.write(str(Parlist[i]))
#			Pfile.close()
			
# ****************************** END WRTITING PARAMETER FILE *************************************
# ****************************** END WRTITING PARAMETER FILE *************************************
# ****************************** END WRTITING PARAMETER FILE *************************************
# ****************************** END WRTITING PARAMETER FILE *************************************
# ****************************** END WRTITING PARAMETER FILE *************************************

# I generate a parameter file from the user imput
		
# ***************************  GENERATING FABADA FILES ******************************************
# ***************************  GENERATING FABADA FILES ******************************************
# ***************************  GENERATING FABADA FILES ******************************************
		a=sqw.getAxis(1)
		qvalues = a.extractValues()		
#		self.log().notice('ispc: '+str(ispc))
#		self.log().notice(str(qvalues))
		
#		self.log().notice(str(qvalues[ispc]))
		
		ispcc = 0
		for ispc in range (SPCmin,SPCmax+1):
			outname="f"+str(ispcc)+".dat"
			outnameV="V"+str(ispcc)+".dat"
			fullpath = os.path.join(directory, outname)
			file = open(fullpath,'w')
			fullpathV = os.path.join(directory, outnameV)
			fileV = open(fullpathV,'w')
#			self.log().notice("Generating: "+str(ispc)+" "+str(ispcc))
#			self.log().notice("Generating: "+str(ispc)+" "+fullpath)
#			self.log().notice("Generating: "+fullpathV)
			ispcc += 1
			
			
			file.write("zvalues 1 "+str(qvalues[ispc])+"\n")
			fileV.write("zvalues 1 "+str(qvalues[ispc])+"\n")
			for idat in range(NB):
				x = sqw.dataX(ispc)[idat]
				y = sqw.dataY(ispc)[idat]
				e = sqw.dataE(ispc)[idat]
				file.write(str(x)+"\t"+str(y)+"\t"+str(e)+"\n")
				x = sqwR.dataX(ispc)[idat]
				y = sqwR.dataY(ispc)[idat]
				e = sqwR.dataE(ispc)[idat]
				fileV.write(str(x)+"\t\t"+str(y)+"\t\t"+str(e)+"\n")
			file.close()
			fileV.close()
# ***************************  END OF GENERATING FABADA FILES ******************************************
# ***************************  END OF GENERATING FABADA FILES ******************************************
# ***************************  END OF GENERATING FABADA FILES ******************************************
# ***************************  END OF GENERATING FABADA FILES ******************************************


# *************************** CALLING FABADA AND EXECUTING IT ************************
# *************************** CALLING FABADA AND EXECUTING IT ************************
		fullpathExe = os.path.join(directory, 'fabada.exe')
		fullpathPfin    = os.path.join(directory, 'Pfin.par')

		subprocess.call([fullpathExe, fullpathPfin])
# *************************** END OF CALLING FABADA AND EXECUTING IT ************************
# *************************** END OF CALLING FABADA AND EXECUTING IT ************************

# *************************** FILLING WORKSPACE WITH RESULTS  ************************
# *************************** FILLING WORKSPACE WITH RESULTS  ************************
# *************************** FILLING WORKSPACE WITH RESULTS  ************************
# OUTPUT FUNCTIONS:
		for ispc in range (NH):
			outname="R"+str(ispc)+".dat"
			fullpath = os.path.join(directory, outname)
#			self.log().notice("Opening: "+fullpath)
			file = open(fullpath,'r')
			line = file.readline()
			i = 0
			while True:
				line = file.readline()
				values = line.split( )
				if line == "":
					break
				
				WSFout.dataX(ispc)[i] = float(values[0])
				WSFout.dataY(ispc)[i] = float(values[2])
				i +=1
			file.close()
		
#Create chi Table
		chiTable = CreateEmptyTableWorkspace(OutputWorkspace="Chi")
		chiTable.addColumn(type="int",name="step") 
		chiTable.addColumn(type="double",name="chi") 
		fullpath = os.path.join(directory, 'chi.dat')
		file = open(fullpath,'r')
		line = file.readline()
		while True:
			line = file.readline()
			values = line.split( )
			if line == "":
				break
				
			chiTable.addRow([ int(values[0]) , float(values[1]) ])
		file.close()
			
#Create Parameter Table
		parameterTable = CreateEmptyTableWorkspace(OutputWorkspace="ParTable")
		parameterTable.addColumn(type="int",name="Number")  # "Detector ID" column required by ApplyCalbration
		parameterTable.addColumn(type="str",name="Name")  # "Detector ID" column required by ApplyCalbration
		parameterTable.addColumn(type="str",name="lg")  # "Detector ID" column required by ApplyCalbration
		parameterTable.addColumn(type="double",name="Parameter")  # "Detector Position" column required by ApplyCalbration
		parameterTable.addColumn(type="double",name="Parjump")  # "Detector Position" column required by ApplyCalbration
		parameterTable.addColumn(type="double",name="Minimum")  # "Detector Position" column required by ApplyCalbration
		parameterTable.addColumn(type="double",name="Maximum")  # "Detector Position" column required by ApplyCalbration
		
		fullpath = os.path.join(directory, 'Pfin.par')
		file = open(fullpath,'r')
		write = False
		while True:
			line = file.readline()
				
			if line == "":
				break
				
			values = line.split( )
			if write == True and len(values) > 3:
#				self.log().notice(str(values[4])+" "+str(values[5])+" "+str(values[6]))
				parameterTable.addRow([ int(values[0]) , str(values[1]) , str(values[2]),float(values[3]) , float(values[4]) , float(values[5]) , float(values[6]) ])
			if "#i" in values:
				write = True
		file.close()
				
#Create MMchain Table
		parameterTable = CreateEmptyTableWorkspace(OutputWorkspace="MMchain")
		parameterTable.addColumn(type="int",name="Steps")  # "Detector ID" column required by ApplyCalbration
		parameterTable.addColumn(type="int",name="DumN")  # "Detector ID" column required by ApplyCalbration
		for i in range(Npar):
			parameterTable.addColumn(type="double",name=Parname[i])  # "Detector ID" column required by ApplyCalbration
		for i in range(Nfunc):
			parameterTable.addColumn(type="double",name="chisq"+str(i))  # "Detector ID" column required by ApplyCalbration
			
		fullpath = os.path.join(directory, 'MMchain.dat')
		file = open(fullpath,'r')
		write = False
		while True:
			roww = []
			line = file.readline()
				
			if line == "":
				break
				
			values = line.split( )
			if write == True:
#				self.log().notice(str(values[4])+" "+str(values[5])+" "+str(values[6]))
				roww.append(int(values[0]))
				roww.append(int(values[1]))
				for i in range(Npar):
					roww.append(float(values[2+i]))
				for i in range(Nfunc):
					roww.append(float(values[2+Npar+i]))
				parameterTable.addRow(roww)
			if "step" in values:
				write = True
		file.close()
			

					

		
# *************************** END OF FILLING WORKSPACE WITH RESULTS  ************************
# *************************** END OF FILLING WORKSPACE WITH RESULTS  ************************
# *************************** END OF FILLING WORKSPACE WITH RESULTS  ************************

AlgorithmFactory.subscribe(Fabada)

