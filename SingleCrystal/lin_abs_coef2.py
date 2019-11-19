#-----------------------------------
#           lin_abs_coef v2.py
#-----------------------------------
from __future__ import print_function

# Program to calculate linear absorption coefficients and density.
# Version 1 requires ISAW. Version 2 is a stand alone script.

# Version 2:
# A.J. Schultz, August 2015

if sys.version_info > (3,0):
    from tkinter import *
    import tkinter.simpledialog as tkSimpleDialog
else:
    from Tkinter import *
    import tkSimpleDialog

import math

class MyDialog(tkSimpleDialog.Dialog):

    def body(self, master):
    
        Message(master, 
            text = (
                "Here are sample inputs of molecular formulas:"
                + "\nFor example, for oxalic acid dihydrate, C2O4H2.2H2O, or C2O6H6, input"
                + "\n     C2 O4 H2 H4 O2    or    C2 O6 H6"
                + "\nFor deuterated oxalic acid dihydrate, input"
                + "\n     C2 O6 D6   or   C2 O6 2H6"
                + "\nFor La2NiO4.2, input"
                + "\n     La2 Ni1 O4.2"
                + "\nFor boron-11 B4C, input"
                + "\n     11B4 C1"
                + "\n\nDescriptions of input:"
                + "\nMolecular formula: The chemical formula input as described above."
                + "\nZ: The number of formula units in the unit cell."
                + "\n   This can be a noninteger value."
                + "\nUnit cell volume: The unit cell volume in units of Angstroms cubed."
                + "\nweight: Crystal weight in milligrams."
                + "\n   If 0, do not calculate crystal volume or radius.\n"),
            # aspect = 200, background = 'yellow').grid(row=0, columnspan=2)
            ).grid(row=0, columnspan=3)
                            
        self.label_text = []
        self.label_text.append( "Molecular formula: " )
        self.label_text.append( "Z: " )
        self.label_text.append( "Unit cell volume: " )
        self.label_text.append( "weight: " )
        self.label_text.append( "Directory for NIST_cross-sections.dat file: " )

        self.entry = []
        for i in range( len( self.label_text) ):
            j = i+1
            Label(master, text = self.label_text[i]).grid(row = j, column = 0, sticky = E)               
            self.entry.append( Entry( master, width = 30) )
            if i == 3: self.entry[i].insert( 0, 0 )
            if i == 4: self.entry[i].insert( 0, '/SNS/software/ISAW/Databases/' )
            self.entry[i].grid( row = j, column = 1, sticky = W )
        
    def apply(self):
    
        self.result = []
        for i in range ( len(self.entry) ):
            self.result.append( self.entry[i].get() )

  
#------------------ Begin ------------------
root = Tk()
root.withdraw()
root.title("lin_abs_coef input")
d = MyDialog(root) 

# Get user input
formulaString = d.result[0]
formulaList = formulaString.split()
numberOfIsotopes = len(formulaList)     # the number of elements or isotopes in the formula        
zParameter = float( d.result[1] )  # number if formulas in the unit cell
unitCellVolume = float( d.result[2] ) # unit cell volume in A^3
# calcRadius = self.getParameter(3).value
calcRadius = True
weight = float( d.result[3] )
XsecDirectory = d.result[4]

sumScatXs = 0.0
sumAbsXs = 0.0
sumAtWt = 0.0

logFileName = 'lin_abs_coef.log'
logFile = open( logFileName, 'w' )
logFile.write('Output from lin_abs_coef.py script:\n\n')

logFile.write('Chemical formula: ' + formulaString + '\n')
logFile.write('Number of formula units in the unit cell (Z): %6.3f\n' % zParameter)
logFile.write('Unit cell volume (A^3): %8.2f\n' % unitCellVolume)

logFile.write('\nCross sections in units of barns ( 1 barn = 1E-24 cm^2)\n')
logFile.write('Absorption cross section for 2200 m/s neutrons (wavelength = 1.8 A)\n')
logFile.write('For further information and references, see ...\ISAW\Databases\NIST_cross-sections.dat\n')

print('\nAtom      ScatXs      AbsXs')	# print headings
print('----      ------      -----')
logFile.write('\nAtom      ScatXs      AbsXs\n')
logFile.write(  '----      ------      -----\n')

# Except for hydrogen, cross-section values are from the NIST web site:
# http://www.ncnr.nist.gov/resources/n-lengths/list.html
# which are from:
# V. F. Sears, Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
# Hydrogen cross-sections are from:
# Howard, J. A. K.; Johnson, O.; Schultz, A. J.; Stringer, A. M.
#	J. Appl. Cryst. 1987, 20, 120-122.
        
filename = XsecDirectory + 'NIST_cross-sections.dat'

# begin loop through each atom in the formula
for i in range(numberOfIsotopes):

    lenAtom = len(formulaList[i])   # length of symbol plus number in formula
    
    # begin test for number of characters in the isotope name or symbol
    for j in range(lenAtom):          
        lenSymbol = lenAtom - j - 1
        if formulaList[i][lenSymbol].isalpha(): break
    lenSymbol = lenSymbol + 1
    
    input = open(filename, 'r')         # this has the effect of rewinding the file
    lineString = input.readline()       # read the first comment line
    while lineString[0] == '#':         # search for the end of the comments block
        lineString = input.readline()

    # Begin to search the table for element/isotope match.
    
    lineList = lineString.split()       # this should be the H atom

    while formulaList[i][0:lenSymbol] != lineList[0]:
        lineString = input.readline()
        lineList = lineString.split()

    scatteringXs = float(lineList[1])   # the total scattering cross section
    absorptionXs = float(lineList[2])   # the true absorption cross section at 1.8 A
    atomicWeight = float(lineList[4])   # atomic weight
    number = float(formulaList[i][lenSymbol:])   # the number of this nuclei in the formula
    
    print('%-5s %10.5f %10.5f' % (lineList[0], scatteringXs, absorptionXs))
    logFile.write('%-5s %10.5f %10.5f\n' % (lineList[0], scatteringXs, absorptionXs))
    
    sumScatXs = sumScatXs + ( number * scatteringXs )
    sumAbsXs = sumAbsXs + ( number * absorptionXs )
    sumAtWt = sumAtWt + ( number * atomicWeight )
    
    input.close()
# end loop

# Calculate the linear absorption coefficients in units of cm^-1
muScat = sumScatXs * zParameter / unitCellVolume
muAbs = sumAbsXs * zParameter / unitCellVolume

# Calculate the density of the crystal in g/cc
density = (sumAtWt / 0.6022) * zParameter / unitCellVolume

# Print the results.
print('\n')
print('The linear absorption coefficent for total scattering is %6.3f cm^-1' % muScat)
print('The linear absorption coefficent for true absorption at 1.8 A is %6.3f cm^-1' % muAbs)
print('The calculated density is %6.3f grams/cm^3' % density)
logFile.write('\n')
logFile.write('The linear absorption coefficent for total scattering is %6.3f cm^-1\n' % muScat)
logFile.write('The linear absorption coefficent for true absorption at 1.8 A is %6.3f cm^-1\n' % muAbs)
logFile.write('\nThe calculated density is %6.3f grams/cm^3\n' % density)

# if calcRadius:
if weight != 0.0:
    crystalVolume = weight / (density)   # sample volume in mm^3
    print('For a weight of %6.3f mg, the crystal volume is %6.3f mm^3' % (weight, crystalVolume))
    logFile.write('\nFor a weight of %6.3f mg, the crystal volume is %6.3f mm^3\n' % (weight, crystalVolume))
    crystalRadius = ( crystalVolume / ((4.0/3.0)*math.pi) )**(1.0/3.0)   # radius in mm
    print('The crystal radius is %6.3f mm, or %6.4f cm' % (crystalRadius, crystalRadius/10.))
    logFile.write('The crystal radius is %6.3f mm, or %6.4f cm\n' % (crystalRadius, crystalRadius/10.))
    # volCalc = (4.0/3.0) * math.pi * crystalRadius**3
    # print 'volCalc = %6.3f' % volCalc

logFile.close()
print('All done!')

