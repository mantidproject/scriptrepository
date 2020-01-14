from __future__ import print_function
import numpy as np
import re, csv
import matplotlib.pyplot as plt

def pretty_print(name, parameter):
  print("===================")
  print("The value of "+name+" is:\n")
  print(parameter)
  print("\n===================\n")

def open_file(filename):
    try:
        thefile = open(filename, 'r')
    except IOError as error:
        print("Problem opening file {file}: {fileerror}".format(file=filename, fileerror=error))
        if thefile != None:
            thefile.close()
            print("Could not open file " + filename + "\n")
    return thefile

def save_dat(filename, data, coldelimiter = "\t"):
    outputfile = open(filename, 'wb')
    line2write = csv.writer(outputfile, delimiter = coldelimiter)
    for line in data:
        line2write.writerow(line)
    del line2write
    outputfile.close()

def save_dict(filename, dictionary, coldelimiter = "\t"):
    outputfile = open(filename, 'wb')
    for key in list(dictionary.keys()):
        line = str(key) + coldelimiter + str(dictionary[key]) + '\n'
        outputfile.write(line)
    del line
    outputfile.close()

def stripheaders(thefile):
    content = []
    for line in thefile:
        m = re.match('[#%]',line)
        if not m:
            content.append(line)
    return content
        
def match_rows_with_numbers(line):
    return re.findall('[0-9Ee.+-]+', line) # careful, this is quite primitive
        # and therfore likely the cause of your problem 
        
def read_dat_file(filename):
    thefile = open_file(filename)
    data = []
    content = stripheaders(thefile)
    for line in content:
        m = match_rows_with_numbers(line) # careful, primitive regex inside 
        if m:
            row = [float(number) for number in m]
            data.append(row)
        m = None
    data = np.array(data)
    return data
   
def make_filename(inputfilename, ending,  separator = '.'):
    '''  '''
    fullname = inputfilename.split(separator,1)
    filename = fullname[0] + ending
   
    return filename

def generate_genx_filenames(inputfilename, howmany = 1):
	fullname = inputfilename.split('.',1)
	filenames = []
	counter = int(fullname[0][-1])
	print(fullname[0][:-1])
	for each in range(howmany):
		counter += 1
		newname = fullname[0][:-1]+str(int(counter))+'.dat'
		filenames.append(newname)
	return filenames
	
	
	
    
def files_to_datadict(filename, list):
    datadict = {}
    for item in list:
	theFile = None
        newFilename = make_filename(filename, item)
        pretty_print("newFilename", newFilename)
        try:
            data = read_dat_file(newFilename)
            datadict[item] = data
        except: pass
      
    return datadict
  
def find_spinfiles(filename, filename2 = None):
    endings = ['.u', '.d', '.uu', '.dd']
    spinstate = files_to_datadict(filename, endings)
    if not spinstate:
        if filename2:
            spinstate['.u'] = read_dat_file(filename) # for those who can't follow useful conventions
            spinstate['.d'] = read_dat_file(filename2)
        else:
            raise Exception('Not enough spin states supplied, could not find two valid files.')
    if '.u' in list(spinstate.keys()):
      spinup = spinstate['.u']
    elif '.uu' in list(spinstate.keys()):
      spinup = spinstate['.uu']
    else: raise Exception('Spinup not found!')
  
    if '.d' in list(spinstate.keys()):
      spindown = spinstate['.d']
    elif '.dd' in list(spinstate.keys()):
      spindown = spinstate['.dd']
    else: raise Exception('Spindown not found!')
    pretty_print("spinup shape", spinup.shape)
    pretty_print("spindown shape", spindown.shape)
    
    return spinup, spindown
    
    
def make_data_finite(data, replacewith = [1e-8,0.1]):
  ''' Check that data and error only contains finite, non-zero numbers and if not replace with default value. Data should have shape (lots, 3), ie x, y, e'''
  cleandata = []
  for line in data:
    for i in [1, 2]:
      testcondition = np.logical_or(np.logical_not(np.isfinite(line[i])), line[i] == 0.0)
      line[i] = np.where(testcondition, replacewith[i-1], line[i])
    cleandata.append(line)
  cleandata = np.array(cleandata)
  return cleandata

def calc_diff(filename1, filename2, newFilename):
  sa1 = read_dat_file(filename1)
  sa2 = read_dat_file(filename2)
  
  difference = sa1 - sa2
  diffData = np.column_stack((sa1[:,0], difference[:,1]))
  
  save_dat(newFilename, diffData)
  
  return diffData
    
    
def sum_spinflip(filenameUD, filenameDU, newFilename):
  '''This sums the two spin-flip data sets. Expects array like structures with shape (rows, 3) ie x,y,e dat files. Saves to newFilename and returns data.'''
  du = read_dat_file(filenameUD)
  ud = read_dat_file(filenameDU)
  pretty_print("the first line ud",ud[0,:])
  pretty_print("the first line du", du[0,:])
  summed = (du+ud)/2
  #recalc the error
  error = np.sqrt(du[:,2]**2+ud[:,2]**2)
  # np.column_stack((x_values, datasummed, errorsummed))
  data = np.column_stack((du[:,0], summed[:,1], error))
  # remove NANs and 0s...
  data = make_data_finite(data, replacewith = [1e-8,0.1])
  # currently have neg numbers below critical edge!
  save_dat(newFilename, data)
  return data
    
def calc_sa(spinup, spindown, factor=1):
  
    theSum = spinup[:,1] + spindown[:,1]
    theDifference = spinup[:,1] - spindown[:,1]
    
    sa = theDifference/theSum * factor
    deltau = spinup[:,2]
    deltad = spindown[:,2]
    error = factor * np.sqrt(deltau**2/theSum**2*(1+theDifference/theSum)**2 + deltad**2/theSum**2*(-1+theDifference/theSum)**2)
    
    sadata = np.column_stack((spinup[:,0], sa, error))
    
    return sadata
    
def make_sa(filename, filename2 = None, factor = 1.0):
   
    spinup, spindown = find_spinfiles(filename, filename2)
    sadata = calc_sa(spinup, spindown, factor)
    
    newFilename = make_filename(filename, '.sa')
    save_dat(newFilename, sadata)
    
    return sadata
  
def genx2sa(filename, filename2 = None):
    
    spinup, spindown = find_spinfiles(filename, filename2)
    sa_data = calc_sa(spinup[:,[0,2,3]], spindown[:,[0,2,3]])
    sa_sim = calc_sa(spinup[:,[0,1,1]], spindown[:,[0,1,1]]) #there is no valid error for the simulation
    sa = np.column_stack((sa_data, sa_sim[:,1]))
    newFilename = make_filename(filename, '.sa')
    save_dat(newFilename, sa)
    
    return sa
    
  
def plotgenxsa(filename, xlim = (0.008, 0.2), ylim = (-1, 1)):
	filename2 = generate_genx_filenames(filename)
	print(filename2)
	sa = genx2sa(filename, filename2[0])
	plt.figure(1,figsize=(25,15))
	plt.ylim(*ylim)
	plt.xlim(*xlim)
	plt.errorbar(sa[:,0], sa[:,1], sa[:,2], fmt='bo', ms=10)
	plt.plot(sa[:,0], sa[:,3], 'r', lw =3)
	plt.ion()
	plt.show()
	

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    list = zip(a,b)
    print(type(list))
    next(b, None)
    return list
	
def theta2Q(theta, lambdas=None):
	'''
	Very simple function that converts a given theta (not 2theta) value 
	to momentum transfer. Theta is in degrees. 
	Use like this:
		theta2Q(0.5) # will print the Q values for a range of wavelengths
		theta2Q(theta=0.5) # as above
		theta2Q(0.7, lambdas=(1.0, 5.0)) #Will print Q for 1.0A and 5.0A
		theta2Q(0.7, 5.0) #Will print Q for theta =0.7deg and lambda=5.0A'''
	if not lambdas:
		lambdas = (1.0, 2.0,4.0,6.0, 8.0, 10.0, 12.0, 14.0)
	try:
		test = lambdas[0] #is it iterable?
	except:
		lambdas = [lambdas] # make it iterable
	print("For theta = "+str(theta))	
	for l in lambdas:
		q = 4*np.pi*np.sin(theta/180*np.pi)/l
		print("lambda: "+str(l)+' -> Q: {:.4f}'.format(q))
			
def Q2theta(Q, lambdas=None):
	'''
	Very simple function that converts a given Q (momentum transfer) 
	to Theta. Theta is in degrees. 
	Use like this:
		Q2theta(0.1) # will print the theta values for a range of wavelengths
		Q2theta(Q=0.1) # as above
		Q2theta(0.06, lambdas=(1.0, 5.0)) #Will print theta for 1.0A and 5.0A
		Q2theta(0.06, 5.0) #Will print theta for Q=0.06A^-1 and lambda=5.0A'''
	if not lambdas:
		lambdas = (1.0, 2.0,4.0,6.0, 8.0, 10.0, 12.0, 14.0)
	try:
		test = lambdas[0] # is it iterable?
	except:
		lambdas = [lambdas] # now it is
	print("For Q = "+str(Q))
	for l in lambdas:
		theta = np.arcsin(l*Q/(4*np.pi))/np.pi*180
		print("lambda: "+str(l)+' -> theta: {:.3f}'.format(theta))

def remove_unwanted_sld_characters(textline):
    # these are very genx sld file specific, the number format is (0.000e+01+0.0000e-01)
    # we need to extract the two numbers
    textline = re.sub('[()]',' ', textline)
    textline = re.sub('e\+','e', textline)
    textline = re.sub('\+', ' ', textline)
    return textline
        
def extract_genx_sld(filename):
    thefile = open_file(filename)
    content = stripheaders(thefile)
    data = [] 
    for line in content:
        line = remove_unwanted_sld_characters(line)
        row = match_rows_with_numbers(line)
        row = [float(number) for number in row]
        data.append(row)
    savewiththisname = make_filename(filename, '_ext.dat',  separator = '.')    
    save_dat(savewiththisname, data)
    return data
        
					