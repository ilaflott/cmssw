#!/bin/sh

for sample in mc data
  do
  for algo in ak
    do
    for sub in Vs Pu
      do
      for radius in 2 3 4 5 6 7
	do
	matchobject="Calo"
	for object in PF Calo
	  do
	  ismc="False"
	  corrlabel="_hiIterativeTracks"
	  if [ $object == "Calo" ]; then
	      corrlabel="_HI"
	  fi
	  if [ $sample == "mc" ]; then
              ismc="True"
          fi

          ismc="False"

	  corrname=`echo ${algo} | sed 's/\(.*\)/\U\1/'`${radius}${object}${corrlabel}
	  genjets="HiGenJets"
          genparticles="hiGenParticles"
	  tracks="hiGeneralTracks"
	  pflow="particleFlowTmp"
	  match=${algo}${sub}${radius}${matchobject}
	  domatch="False"
#	  matchobject="PF"

	  cat templateSequence_cff.py.txt \
	      | sed "s/ALGO_/$algo/g" \
	      | sed "s/SUB_/$sub/g" \
	      | sed "s/RADIUS_/$radius/g" \
	      | sed "s/OBJECT_/$object/g" \
	      | sed "s/SAMPLE_/$sample/g" \
	      | sed "s/CORRNAME_/$corrname/g" \
	      | sed "s/MATCHED_/$match/g" \
	      | sed "s/ISMC/$ismc/g" \
              | sed "s/GENJETS/$genjets/g" \
              | sed "s/GENPARTICLES/$genparticles/g" \
              | sed "s/TRACKS/$tracks/g" \
              | sed "s/PARTICLEFLOW/$pflow/g" \
	      | sed "s/DOMATCH/$domatch/g" \
	      > $algo$sub$radius${object}JetSequence_${sample}_cff.py	  	  
	done
      done
    done
  done
done





