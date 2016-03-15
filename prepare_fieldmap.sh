#!/bin/sh

# usage: prepare_fieldmap.sh bids-dir bids-sub delta-TE

# check if subject folder exists 
if [ ! -d "$1/derivatives/$2" ]; then
	mkdir $1/derivatives/$2
fi

# check if fmap folder exists
if [ ! -d "$1/derivatives/$2/fmap" ]; then
	mkdir $1/derivatives/$2/fmap
fi

# extract second echo from magnitude image
$FSLDIR/bin/fslroi $1/$2/fmap/$2_magnitude $1/derivatives/$2/fmap/$2_magnitude_e2 1 1

# bet magnitude image and erode
$FSLDIR/bin/bet $1/derivatives/$2/fmap/$2_magnitude_e2 $1/derivatives/$2/fmap/$2_magnitude_e2_brain -f 0.5 -g 0
$FSLDIR/bin/fslmaths $1/derivatives/$2/fmap/$2_magnitude_e2_brain -ero $1/derivatives/$2/fmap/$2_magnitude_e2_brain_ero

# prepare field map
$FSLDIR/bin/fsl_prepare_fieldmap SIEMENS $1/$2/fmap/$2_phasediff $1/derivatives/$2/fmap/$2_magnitude_e2_brain_ero $1/derivatives/$2/fmap/$2_fieldmap $3
