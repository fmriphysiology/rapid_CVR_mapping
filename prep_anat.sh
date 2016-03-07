#!/bin/sh

# usage: prep_anat.sh bids-dir bids-sub 

# check if subject folder exists 
if [ ! -d "$1/derivatives/$2" ]; then
	mkdir $1/derivatives/$2
fi

# check if fmap folder exists
if [ ! -d "$1/derivatives/$2/anat" ]; then
	mkdir $1/derivatives/$2/anat
fi

# bet anatomical
$FSLDIR/bin/bet $1/$2/anat/$2_T1w $1/derivatives/$2/anat/$2_T1w_brain -f 0.5 -g 0
