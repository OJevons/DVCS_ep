#!/bin/bash
###########################################################
# Create DVCS filelist from S3 folder                     #
###########################################################

###########################################################
# Help message                                            #
###########################################################
Help()
{
    # Display help
    echo "Create list of ePIC simulation campaign files by campaign date and beam settings."
    echo "NOTE: This script MUST be run within eic-shell to work."
    echo
    echo "Syntax:  getFilelist.sh [-b|c|e|h]"
    echo "options:"
    echo "b     Beam setting (hiAcc, hiDiv)."
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

#echo "Beam setting: $BeamSet"
#echo "Beam energy: $Energy"
#echo "Campaign: $Camp"

# Check for valid inputs: beam setting
if [ $BeamSet == "X" ]
then
    echo "Must declare beam setting."
    exit
elif [ $BeamSet != "hiDiv" ] && [ $BeamSet != "hiAcc" ]
then
    echo "Invalid beam setting."
    exit
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


mc ls S3/eictest/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/$Energy | grep $BeamSet | sed "s,.*DVCS.,root://dtn-eic.jlab.org//work/eic2/EPIC/RECO/$Camp/epic_craterlake/EXCLUSIVE/DVCS_ABCONV/$Energy/DVCS.," > ../DVCS_ep/filelists/inputFileList_ePIC_"$Camp"_"$Energy"_"$BeamSet".list

