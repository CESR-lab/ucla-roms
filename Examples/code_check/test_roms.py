# Python script that takes the log file of a benchmark run and compares it
# to the latest result log for that example.
# -> compares the diagnostic numbers to ensure they haven't changed.
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

# -----------------------------------------------------------------

# ******* USER INPUTS START ************
# Tests & file names of roms terminal output logs
test_names = ['WEC_real','Rivers_ana'] # ,'Rivers_ana'
ntests = np.size(test_names)
# -- WEC
filename_BM  = os.path.join(os.getcwd(),'benchmarks/wec/sample_wec_bm.log')
filename_res = os.path.join(os.getcwd(),'benchmarks/wec/sample_wec_test.log')
filenames = [ filename_BM, filename_res ]
#print(filenames)
# -- Rivers_ana
filename_BM  = os.path.join(os.getcwd(),'benchmarks/river/river_ana_bm.log')
filename_res = os.path.join(os.getcwd(),'benchmarks/river/river_ana_test.log')
filenames = np.append(filenames,[ filename_BM, filename_res ])
filenames = np.reshape(filenames,(ntests,2))
#print(np.shape(filenames))
# -- netcdf test of his, avg & rst
fname_ncdf_BM  = os.path.join(os.getcwd(),'benchmarks/netcdf/wec_netcdf_bm.log')
fname_ncdf_res = os.path.join(os.getcwd(),'benchmarks/netcdf/wec_netcdf_test.log')
fnames_ncdf = [ fname_ncdf_BM, fname_ncdf_res ]
# ******* USER INPUTS END **************

# --------------------------------------------
# Diagnostic terms
fields = ['KINETIC_ENRG','BAROTR_KE','MAX_ADV_CFL','MAX_VERT_CFL']
diags = np.zeros((ntests,2,4)) # 2 for bm vs new, 4 for fields
diffs = np.zeros((ntests,4)) # Differences between bm vs new fiels - 4 for fields

# Generate log file for diagnostic checks from this script
file_diag = open("code_check.log","w")
file_diag.write("ROMS code check log:\n")

# Loop through different test cases:
for test in range(ntests):

	text='Checking test: '+test_names[test]
	print(text)
	file_diag.write('\n'+text) # write("\n"+text)
	# -- read number of timesteps
	line_num = search_string_in_file(filenames[test,0],'ntimes')
	line = open(filenames[test,0], 'r').readlines()[line_num[0]-1] # -1 for 0 counting
	nstps = int(line[18:23]) + 1 # use +1 as diagnostics start from step 0
	#print('nstps='+nstps)
	# -- find line where diagnostics start
	lstart = search_string_in_file(filenames[test,0],'STEP')

	# Loop through BM values and new values
	for m in range(2):

		# Loop through time steps & sum diagnostics
		for t in range(nstps):

			# read timestep line
			line = open(filenames[test,m], 'r').readlines()[lstart[0]+t]
			# read KINETIC_ENRG
			val = line[18:34]
			# For some reason that I could not solve some, the .log numbers in scientific 
			# notation are missing the 'E' exponent, and hence python does not recognize 
			# the number. Numbers are written in Diag.F. Hence need to add the 'E'.
			diags[test,m,0] += float(val[:13] + "E" + val[13:])
			#print('diags[0]',diags[0])    
			# read BAROTR_KE
			val = line[35:50]
			diags[test,m,1] += float(val[:12] + "E" + val[12:])
			#print('diags[1]',diags[1])
			# read MAX_ADV_CFL
			diags[test,m,2] += float(line[52:66])
			#print('diags[2]',diags[2])
			# read MAX_VERT_CFL
			diags[test,m,3] += float(line[68:82])
			#print('diags[3]',diags[3])
			
			# Sanity check to confirm correct read of values as totals only
			# are hard to spot check
			if t==1 and m==1:
				file_diag.write('\n'+'Diagnostics for 1st timestep of '+str(nstps-1)+' total steps:\n') # -1 as +1 before
				file_diag.write( str(t)+' '+str(diags[test,m,:])+'\n' )				
	
	# Confirm results
	file_diag.write('RESULTS for: '+test_names[test]+'\n')
	for t in range(4):
		file_diag.write(fields[t]+':\n')
		file_diag.write(str(diags[test,0,t])+' - BM\n')  # debug ' = diags[0,'+str(t)+'] = ', 
		file_diag.write(str(diags[test,1,t])+' - new\n') # debug ' = diags[1,'+str(t)+'] = ',
		diffs[test,t] = diags[test,1,t]-diags[test,0,t]	
		file_diag.write(str(diffs[test,t])+' = difference\n')

	# Confirm outcome to user:
	if np.all(diffs[test,:]==0.0):
		text=test_names[test]+' test is correct! --> diagnostics match the benchmark.\n'
		print(text)
		file_diag.write(text)
	else:
		text='ERROR! --> test does not match benchmark result.\n'
		print(text)
		file_diag.write(text)

# ----------------------------------------------------------
# Confirm netcdf functionality working for 3 types:
# his (history) , avg (averages), and rst (restart) files
#
# Crudely checks that run has outputted the desired number
# of steps. Test that results are correct is covered in
# tests above. Didn't include the netcdf output above because
# the printing of 'wrt_his: ...' messes up the loop of counting
# the diagnostic values.

# netcdf terms
types = ['history','restart','averages']
outputs   = np.zeros((2,3)) # 2 for bm vs new, 3 for types
out_diffs = np.zeros((3)) # Differences between bm vs new fiels - 3 for types

text='Checking netcdf: outputting of history, restart & averages'
print(text)
file_diag.write('\n'+text+'\n')

# Loop through BM values and new values
for m in range(2):

	# read number of timesteps
	line_num = search_string_in_file(fnames_ncdf[m],'Records written:')
	line = open(fnames_ncdf[m], 'r').readlines()[line_num[0]-1] # -1 for 0 counting
	#print('line= ' + line)
	outputs[m,0] = int(line[42:44]) # Number of history steps written
	#print('history= '  + str(outputs[m,0]))
	outputs[m,1] = int(line[56:58]) # Number of restart steps written
	#print('restart= '  + str(outputs[m,1]))
	outputs[m,2] = int(line[71:73]) # Number of history steps written
	#print('averages= ' + str(outputs[m,2]))

out_diffs = outputs[0,:]-outputs[1,:]

# Confirm outcome to user:
if np.all(out_diffs==0.0):
	text='Netcdf output works! --> numbers of history, restart & averages records match the benchmark.\n'
	print(text)
	file_diag.write(text)
else:
	text='ERROR! --> Netcdf output not working.\n'
	print(text)
	file_diag.write(text)	
	
	
# ----------------------------------------------------------
# Confirm overall outcome to user:
if np.all(diffs==0.0) and np.all(out_diffs==0.0):
    text='\nCODE CORRECT! --> All diagnostics & netcdf outputting match the benchmarks. \nSee code_check.log file for details.'
    print(text)
    file_diag.write(text)
else:
    text='\nERROR! --> Tests do not match benchmark results. \nEnsure using correct benchmark for machine. Deafult machine is Maya. \nSee code_check.log file for details.'
    print(text)
    file_diag.write(text)


file_diag.close() # close code check log file



