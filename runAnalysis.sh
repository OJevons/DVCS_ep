#!/bin/bash
###########################################################
# Run full DVCS analysis from command line                #
###########################################################

###########################################################
# Help message                                            #
###########################################################
Help()
{
    # Display help
    echo "Run ePIC ep DVCS analysis on simulation campaign files."
    echo "NOTE: This script MUST be run within eic-shell to work."
    echo
    echo "Syntax:  runAnalysis.sh [-b|c|e|h]"
    echo "options:"
    echo "b     Beam setting (hiAcc, hiDiv, TEST)."
    echo "c     Campaign to look at (23.08.0, 23.10.1, etc.)."
    echo "e     Beam energy configuration (\"5x41\", \"10x100\", \"18x275\")."
    echo
    echo "Help message is run if no options are provided."
}

###########################################################
# Main script                                             #
###########################################################

# Print help message if no arguments are provided
if [[ -z $1 ]];
then
    Help
    exit
fi

# Set initial variables
BeamSet="X"
Camp="X"
Energy="X"


# Parse over all options
while getopts h:b:c:e: option; do
    case $option in
	h) # Print help message
	    Help
	    exit;;
	b) # Set beam setting
	    BeamSet=$OPTARG;;
	c) # Set campaign to use
	    Camp=$OPTARG;;
	e) # Set beam energy
	    Energy=$OPTARG;;
       \?) # Invalid options
	    echo "Error: Invalid flag"
	    exit;;
    esac
done


# Check for input: beam setting
# If no input, then set $BeamSet to "TEST"
# Don't worry about invalid inputs here. Analysis code defaults to a test case
if [ $BeamSet == "X" ]
then
    $BeamSet = "TEST"
fi
# Check for input: campaign
if [ $Camp == 'X' ]
then
    echo "Must declare campaign to use."
    exit
fi
# Check for valid inputs: beam energy
if [ $Energy == "X" ]
then
    echo "Must declare beam energy."
    exit
elif [ $Energy != "5x41" ] && [ $Energy != "10x100" ] && [ $Energy != "18x275" ]
then
    echo "Invalid beam energy."
    exit
fi

# Run analysis
root -q 'run_ePIC_DVCS.C("'$Camp'","'$Energy'","'$BeamSet'")'
# Run plotting macro
root -q 'DVCSPlots.C("'$Camp'","'$Energy'","'$BeamSet'")'
