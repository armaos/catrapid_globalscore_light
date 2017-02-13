#!/bin/bash
# USAGE: bash parallelize.sh <number of jobs> <protein library file> <rna library file>

NUMJOBS=$1

NUMPROTS=`wc -l $2 | awk '{print $1}'`;

# Protein file splitting
grep -v "#" $2 | awk '{sp=(int(("'$NUMPROTS'"/10)/"'$NUMJOBS'"+1)*10);a[NR]=$0}END{f=1; filename=sprintf("list_prot.%d", f);c=1;for(i=1;i<=NR;i++){if(c<=sp){print a[i] >> filename;c++} else{c=1; f++;filename=sprintf("list_prot.%d", f); print a[i] >> filename;c++} }}'

for((i=1;i<=NUMJOBS;i++)); do

	if [[ -s list_prot.$i ]]; then 
		cp -r template dir$i
		mv list_prot.$i dir$i/prot/
		cp $3 dir$i/rna/
		cd dir$i
		bash runlibrary.both.sh prot/list_prot.$i $3 $i &
		cd ..
	fi

done

wait

if [[ -s "pre-compiled/out.merged.txt" ]]; then
	rm pre-compiled/out.merged.txt
fi

if [[ -s "protein.errors.log" ]]; then
	rm protein.errors.log
fi

if [[ -s "rna.errors.log" ]]; then
	rm rna.errors.log
fi

for((i=1;i<=NUMJOBS;i++)); do
	if [[ -s dir$i/pre-compiled/out.$i.txt ]]; then
		awk '{printf "%s %s\t\t%.2f\t%.2f\n", $1, $2, $3, $4}' dir$i/pre-compiled/out.$i.txt >> pre-compiled/out.merged.txt
		cat dir$i/protein.errors >> protein.errors.log
		cat dir$i/rna.errors >> rna.errors.log
		rm -fr dir$i
	fi
done
