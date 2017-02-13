#!/bin/bash

echo "Protein sequence processing and fragmentation"
date +"%m-%d-%y %r"
protein_name=$(echo ${1##*/})
rna_name=$(echo ${2##*/})
./flip -u $1
sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $1 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){printf "%s", toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, toupper($2)}' >  ./protein/outfile

pr=`cat ./protein/outfile | awk '(NR==1){s=length($2)/50; if(s<=25){s=25}   if(s>=375){s=375} print int(s)}'`
cd catrapid-fragmentation
if [ ! -s tmp ]; then
  mkdir tmp/
fi
bash bin/cutter.sh ../protein/outfile   $pr >  ../protein/outfile.fr
rm tmp/*
cd ..

echo "Protein library generation"
date +"%m-%d-%y %r"

#pythonic way
cd library_generator_python
  sed 's/|/_/g' ../protein/outfile.fr > protein.fr.txt
  python library_generator.pyc protein.fr.txt "protein"
  python library_checker.pyc outs
  mv outs/protein.fr.txt.protein.lib ../interactions.U/combine_parallel/prot/protein.lib
  rm -fr ./tmp_protein
cd ..


echo "RNA sequence processing and fragmentation"
date +"%m-%d-%y %r"

./flip -u $2
sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $2 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){gsub(/[Uu]/, "T", $1); printf "%s",toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, $2}' >  ./rna/outfile2

nt=`cat ./rna/outfile2  | awk '(NR==1){s=length($2)/50; if(s<=25){s=25}   if(s>=750){s=750} print int(s)}'`
cd catrapid-fragmentation
if [ ! -s tmp ]; then
  mkdir tmp/
fi
bash bin/cutter.sh ../rna/outfile2  $nt >  ../rna/outfile2.fr
rm tmp/*
cd ..

cd library_generator_python
  sed 's/|/_/g' ../rna/outfile2.fr > rna.fr.txt
  python library_generator.pyc rna.fr.txt "rna"
  python library_checker.pyc outs
  mv outs/rna.fr.txt.rna.lib ../interactions.U/combine_parallel/rna/rna.lib
  rm -fr *.ss
  rm -fr ./tmp_rna
cd ..

cd interactions.U/
    cd combine_parallel/
      if [ ! -s pre-compiled ]; then
        mkdir pre-compiled/
      fi
      python multiplier.pyc 10 "prot" "rna"
cd ../../

cd filter
echo "Global Score computing"
date +"%m-%d-%y %r"

awk '{print "protein_"$1,"rna_"$2,$3,$4}' ../interactions.U/combine_parallel/pre-compiled/* > interactions.txt
bash start.sh interactions.txt -1 > processed.txt
awk '{printf "%.2f\n", ($1+1)/2}' processed.txt > ../outputs/$protein_name.$rna_name.GS.txt

cd ..
