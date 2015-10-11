#!/usr/bin/python

def create_submit_all_script():
    """Do an os walk to find all SLURM scripts that should be submitted."""
    import os
    import sys
    #
    def run_bash(cmd):
        # Runs a bash command and returns the output
        import subprocess
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out = p.stdout.read().strip()
        return out  # This is the stdout from the shell command
    #
    absolute_path = os.path.join(os.getcwd())
    #
    all_slurm_files = []
    all_slurm_file_folders = []
    # Collect all SLURM files to one big list
    for (dirpath,dirnames,filenames) in os.walk(absolute_path):
        for filename in filenames:
            if filename.endswith('.cmd'):
                all_slurm_files.append(filename)
                all_slurm_file_folders.append(os.path.join(absolute_path,dirpath))
    #
    #for i in range(len(all_slurm_files)):
    #    print(all_slurm_file_folders[i],all_slurm_files[i])
    #
    # Finally, create a bash script for submitting all of the calculations:
    #
    script_name = os.path.join(absolute_path,'submit_all')
    submit_script = open(script_name,"w")
    #for_loop_line = 'for i in '
    #for i in range(len(all_slurm_files)):
    #    for_loop_line += '{0} '.format(all_slurm_files[i])
    submit_lines  = '#!/bin/bash'+'\n'
    submit_lines += '\n'
    #submit_lines += for_loop_line + '\n'
    #submit_lines += 'do' + '\n'
    #submit_lines += 'cd $i; sbatch slurm.sh;cd ..' + '\n'
    #submit_lines += 'done' + '\n'
    #
    for i in range(len(all_slurm_files)):
        submit_lines += 'cd {0}; sbatch {1}\n'.format(all_slurm_file_folders[i],all_slurm_files[i])
    #
    for line in submit_lines:
        submit_script.write(line)
    submit_script.close()
    # Give the script executability privileges
    cmd_line = 'chmod +x {0}'.format(script_name)
    run_bash(cmd_line)
    #
    #print(len(all_slurm_files))
    return

if __name__ == "__main__":
    create_submit_all_script()