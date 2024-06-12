# a script to produce events with pythia cards produced by produce_cards.py
# performs a scan of lifetimes and masses for the higgs and vector portal models

import subprocess

def runCommand(command):
	# Start the subprocess
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	
	# Read the output line by line
	while True:
	    output = process.stdout.readline()
	    if output == '' and process.poll() is not None:
	        break
	    if output:
	        print(output.strip())
	    
	# Capture any remaining output after the process has finished
	remaining_output, errors = process.communicate()
	
	if remaining_output:
	    print(remaining_output.strip())
	
	if errors:
	    print("Errors:")
	    print(errors.strip())


# set your card path
card_path="/Users/karri/Documents/ATLAS/darkQCD/dark_showers_tool/cards"
hepmc_path="/Users/karri/Documents/ATLAS/darkQCD/dark_showers_tool/hepmc"
pythia8="/opt/homebrew/Cellar/pythia/8.310/share/Pythia8"

portal="higgsportal"
production="wh"
xio = 1.0
xil = 1.0
masses = [1,2,3,5,10]
ctaus = [3,10,30,100] # mm
production = "wh" # ggf or zh

for mass in masses: 
	for ctau in ctaus: 

	sample=f"{production}_{portal}_m-{mass}_xio-{xio}_xil-{xil}_ctau-{ctau}"
	
	# to do - if not then stop
	# check that pythia example main42 is compiled
	# check that input hepmc file exists
	# check that output directory exists
	
	command = f"{pythia8}/examples/main42 {card_path}/{sample}.cmnd {hepmc_path}/{sample}.hepmc"

	runCommand(command)

