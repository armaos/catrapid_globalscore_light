#!/bin/bash
mkdir ./tmp
mkdir ./results
mkdir ./data
mkdir ./database

#lenP=`echo $1 | sed 's/\.prot\.lib//g;s/rand\.//g' | awk '{print $1}'`
#lenR=`echo $2 | sed 's/rand\.//g;s/\.rna\.lib//g'  | awk '{print $1}'`

touch protein.errors rna.errors

#Protein
for i1 in `cat $1 | grep -v "#" | awk '{print $1}' | uniq | head -200000`; do 

	grep -w $i1 $1  | grep -v "#" | awk '{print $1}' | uniq > protein.list.txt 
	grep -w $i1 $1  | grep -v "#" | awk '{for(i=2;i<=51;i++){printf "\t%4.2f", $i} printf "\n"}' > data/bait.posi

	if [[ `wc -l data/bait.posi | awk '{print $1}'` -eq 10 ]]; then
		# RNA
		for i2 in `cat $2 | grep -v "#" | awk '{print $1}' | uniq | head -100000`; do 

			grep -w $i2 $2 | grep -v "#" | awk '{for(i=2;i<=51;i++){printf "\t%4.2f", $i} printf "\n"}' > data/prey.posi
			
			if [[ `wc -l data/prey.posi | awk '{print $1}'` -eq 10 ]]; then

				bash ./bin/start.modified.sh $i2
				
			else
			
				echo "$i2" >> rna.errors
				
			fi

		done >> pre-compiled/out.$3.txt
		
	else
	
		echo "$i1" >> protein.errors
	
	fi

done

#rm positives.txt negatives.txt protein.list.txt
#rm -fr database
#rm -fr data
#rm -fr results
#rm -fr tmp
