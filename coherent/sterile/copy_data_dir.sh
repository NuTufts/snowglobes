#!/bin/bash

INDIR=$1
OUTDIR=$2

cp -r $INDIR/channels $OUTDIR/
cp -r $INDIR/fluxes   $OUTDIR/
cp -r $INDIR/glb      $OUTDIR/
cp -r $INDIR/xscns    $OUTDIR/
cp -r $INDIR/effic    $OUTDIR/
cp -r $INDIR/smear    $OUTDIR/
cp $INDIR/detector_configurations.dat $OUTDIR/
