# Python script that takes the log file of a benchmark run and compares it
# to the newly run result log for that example.
# -> compares the sum of diagnostic numbers to ensure they haven't changed.
# -> roms run uses own diag.F to get full precision output to compare.
# -> only useful when you expect result to be exactly the same.
#
# Add the appropriate .log file in the USER INPUTS section.
#
# Old roms code did a similar test using 'etalon' as can still be seen
# in diag.F.

################################

# Imports:
import os
import sys
import numpy as np


# Define function to use:
def search_string_in_file(file_name, string_to_search):
    line_number = 0
    list_of_results = []
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                list_of_results.append((line_number))                
    # Return list of tuples containing line numbers
    return list_of_results
    
def search_string_in_file_yesno(file_name, string_to_search):
    line_number = 0
    found_string = False
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            line_number += 1
            if string_to_search in line:
                # If yes, then add the line number & line as a tuple in the list
                found_string=True
                break                
    
    return found_string

# -----------------------------------------------------------------

# ******* USER INPUTS START ************
# Tests & file names of roms terminal output logs
ntests = 1 # np.size(test_names)
# -- WEC
#filename_BM  = os.path.join(os.getcwd(),'benchmark.result')
filename_BM  = os.path.join(os.getcwd(),str(sys.argv[1]))  # read in with argument (maya or expanse)
filename_res = os.path.join(os.getcwd(),'test.log')
filenames = [ filename_BM, filename_res ]
#print(filenames)
# ******* USER INPUTS END **************

# --------------------------------------------
# Diagnostic terms
fields = ['KINETIC_ENRG','BAROTR_KE','MAX_ADV_CFL','MAX_VERT_CFL']
diags = np.zeros((ntests,2,4)) # 2 for benchmark vs new, 4 for fields
diffs = np.zeros((ntests,4)) # Differences between benchmark vs new fields - 4 for fields

# Generate log file for diagnostic checks from this script
file_diag = open("code_check.log","w")
file_diag.write("ROMS code check log:\n")

# Loop through different test cases:
for test in range(ntests):

	text='  checking diagnostic values: ' # +test_names[test]
	print(text)
	file_diag.write('\n'+text) # write("\n"+text)
	# -- read number of timesteps
	line_num = search_string_in_file(filenames[0],'ntimes')   # find line number of ntimes in text
	line = open(filenames[0], 'r').readlines()[line_num[0]-1] # -1 for 0 counting
	nstps = int(line[18:23]) + 1 # use +1 as diagnostics start from step 0
	#print('nstps=',nstps)

	# Loop through Benchmark log values and new log values
	for m in range(2):
	
	  # -- find line where diagnostics start (within loop in-case terminal output changes between commits)
	  lstart = search_string_in_file(filenames[m],'STEP')
	  #print('lstart=',lstart)   # debug
	
	  # Loop through time steps & sum diagnostics
	  iline = lstart[0]
	  dline = 0
	  while dline < nstps :
	    line = open(filenames[m],'r').readlines()[iline]
	    if ((len(line) > 94) and (len(line.split())==5)) :
	      diags[test,m,0] += float(line[ 4:26])     # read Kinetic Energy
	      diags[test,m,1] += float(line[27:49])     # read barotropic KE
	      diags[test,m,2] += float(line[50:72])     # read MAX_ADV_CFL        
	      diags[test,m,3] += float(line[73:95])     # read MAX_VERT_CFL
	      dline += 1
	    iline += 1

#	for t in range(nstps):
#		#print('t=',t)   # debug
#		# read timestep line
#		line = open(filenames[m], 'r').readlines()[lstart[0]+t]  # need the [0] even though scalar
#		  
#		diags[test,m,0] += float(line[ 4:26])          # python's float already double precision (64 bit)
#		#print('diags[test,m,0]',diags[test,m,0])      # read KINETIC_ENRG
#		diags[test,m,1] += float(line[27:49])          # read BAROTR_KE
#		#print('diags[test,m,1]',diags[test,m,1])      
#		diags[test,m,2] += float(line[50:72])          # read MAX_ADV_CFL        
#		#print('diags[test,m,2]',diags[test,m,2])
#		diags[test,m,3] += float(line[73:95])          # read MAX_VERT_CFL
#		#print('diags[test,m,3]',diags[test,m,3])
#		
#		# Sanity check to confirm correct read of values as totals only
#		# are hard to spot check
#		if t==1 and m==1:
#			file_diag.write('\n'+'Diagnostics for 1st timestep of '+str(nstps-1)+' total steps:\n') # -1 as +1 before
#			file_diag.write( str(t)+' '+str(diags[test,m,:])+'\n' )				
	
	# Confirm results
	file_diag.write('RESULTS: \n')
	for t in range(4):
		file_diag.write(fields[t]+':\n')
		file_diag.write(str(diags[test,0,t])+' - BM\n')  # debug ' = diags[0,'+str(t)+'] = ', 
		file_diag.write(str(diags[test,1,t])+' - new\n') # debug ' = diags[1,'+str(t)+'] = ',
		diffs[test,t] = diags[test,1,t]-diags[test,0,t]	
		file_diag.write(str(diffs[test,t])+' = difference\n')

	# Confirm outcome to user:
	if np.all(diffs[test,:]==0.0):
		text='    test is correct! --> diagnostics match the benchmark.'
		print(text)
		file_diag.write(text)
	else:
		text='    ERROR! --> test does not match benchmark result.\n'
		print(text)
		file_diag.write(text)

# ----------------------------------------------------------
# Confirm netcdf functionality:
# - crudely checks output was successful

# netcdf terms
types = ['history','restart','averages']
outputs   = np.zeros((2,3)) # 2 for bm vs new, 3 for types
ncdf_diffs = np.zeros((3)) # Differences between bm vs new fiels - 3 for types

text='  checking netCDF output:'
print(text)
file_diag.write('\n'+text+'\n')

# only need to check in new result if wrote history was successful
found_str = False
found_str = search_string_in_file_yesno(filenames[1],'wrote history')

# Loop through BM values and new values
#for m in range(2):

	# read number of timesteps
	#found_str[m] = search_string_in_file_yesno(filenames[m],'wrote history') # supressed write history so only write it once
	#line = open(filenames[m], 'r').readlines()[line_num[0]-1] # -1 for 0 counting
	#print('line= ' + line)
	#outputs[m,0] = int(line[42:44]) # Number of history steps written
	#print('history= '  + str(outputs[m,0]))
	#outputs[m,1] = int(line[56:58]) # Number of restart steps written
	#print('restart= '  + str(outputs[m,1]))
	#outputs[m,2] = int(line[71:73]) # Number of history steps written
	#print('averages= ' + str(outputs[m,2]))

#ncdf_diffs = outputs[0,:]-outputs[1,:]

# Confirm outcome to user:
if (found_str):
#	text='    netCDF output works! --> numbers of history records matches the benchmark.'
	text='    netCDF output works! --> wrote history string found at end of run.'
	print(text)
	file_diag.write(text)
else:
	text='    ERROR! --> Netcdf output not working.\n'
	print(text)
	file_diag.write(text)	
	
	
# ----------------------------------------------------------
# Confirm overall outcome to user:
if np.all(diffs==0.0) and (found_str):
    text='  CODE CORRECT! --> All diagnostics & netcdf outputting match the benchmarks. \n  see code_check.log file for details.\n'
    print(text)
    file_diag.write(text)
    retval=0                                     # let shell script know check was successful
else:
    text='\n  ERROR! --> Tests do not match benchmark results.\n  See code_check.log file for details.'
    print(text)
    file_diag.write(text)
    retval=1                                     # let shell script know check failed


file_diag.close() # close code check log file
sys.exit(retval)



