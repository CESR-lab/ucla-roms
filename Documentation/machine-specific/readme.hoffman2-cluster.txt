Additional instructions for ROMS to compile on Hoffman2 cluster at UCLA IDRE.

1) Copy the bashrc file stored on the repo as your new .bashrc (save your old one somewhere).
   Edit the line for your ROMS_ROOT. (don't forget to 'source ~/.bashrc' afterwards):
   Documentation/machine-specific/bashrc-cshrc-files/bashrc-hoffman2-202403

2) Edit NHMG/src/Makefile - uncomment line: FC = mpiifort.

3) Uncomment the two mpiifort lines in src/Makedefs.inc.

4) Follow all the usual compile steps in the readme's.

General Hoffman2 information:
  Hoffman2 is a UCLA IDRE cluster. You need to be part of a group to access and run on Hoffman2. 

How to run ROMS on Hoffman2:
  Sample run scripts
  Documentation/machine-specific/Job_submission_scripts
  Submit job: qsub <your_run_script.cmd> (submit.cmd in the example)
  Job status: qstat -u <your_username>
  Cancel job: qdel <job_number> or qdel -u <your_username> (to delete all your jobs)

Running an example:
  Compile Examples/Pipes_ana/
  Copy the Hoffman2 scripts from Documentation/machine-specific/Job_submission_scripts to Examples/Pipes_ana/
  Edit submit.cmd and do_roms_hoffman.bash to your own directories
  Enter the command: qsub submit.cmd
  Check the roms.joblog to see if it succeeded

Hoffman2 job submission options:
  Every job must be submitted with requested number of cores, memory, and time. 
  #  Resources requested
  #$ -pe dc* 6                # requested number of nodes, 6 here
  #$-l h_data=4G              # h_data is memory requested, 4 GB here
  #$-l h_rt=1:00:00           # h_rt is time requested, 1 hour here
  #$-l highp                  # highp keyword is to use only your group's owned nodes; remove if using public nodes
  #$-l arch=intel-gold-6140   # arch is to specify a specific node architecture; remove if unknown or unsure
			      # only useful if you want to use certain nodes or run many jobs simultaneously

  Putting it all together would look like this in the job submission
  #$ -pe dc* 6
  #$-l h_data=4G,h_rt=28:00:00,arch=intel-gold-6140,highp

Hardware:
  Hoffman2 has shared/public nodes you can use for up to 24 hours. You can use your group's owned nodes for unlimited time and queue
multiple jobs in a row on those nodes. If your group's nodes are not being used, they will be lent out to other group as public nodes
for use up to 24 hours. If trying to use your group's owned nodes and your job doesn't start immediately, they are likley being
lent out and you will have to wait until they are freed. 

Transferring data:
  If you have lots of large files to transfer between your personal machine / campus server,
  consider using this web portal: globus.org.  You log in with your XSEDE/ACCESS or Hoffman2 account, then it provides a useful
  webpage GUI to transfer all your files.
  In the GUI use:
  Collection = Official UCLA Hoffman2 Cluster
  home path = /u/home/$first_letter_of_username$/$username$ (to get to your expanse home directory)
  scratch path = /u/scratch/$first_letter_of_username$/$username$ 
  rsync can be used but is generally much slower.

Resources and documentation:
https://www.hoffman2.idre.ucla.edu/
