# runs like : sh assessor.sh c.u.txt 41 -1
# where 41 is the model and +1 / -1 is weighted / uniform

# if the command line contains +1 or -1 it formats the set, otherwise it uses a different routine with letters U and W

# marks a file that will be anyway deleted
echo "del" > input.txt



s=`echo $3  | awk 'BEGIN{s=0} ($1~/1/){s=1} END{print s}'`

# preformatted
if [ $s == 0 ]; then
echo "# pre-formatted set"
for j in `awk '($NF~/\'$4'/){print $1}' $1 |  sed 's/NM_/NM/g' |  sed 's/NP_/NP/g' | sed 's/XM_/XM/g'  | sed 's/_/ /g' | awk '{print $1}' | sort -u | sed 's/NM/NM_/g' | sed 's/NP/NP_/g'| sed 's/XM/XM_/g'`;
do 
grep $j $1 | sort -k 2 |  uniq | tail -1 > in.tmp;  
sh start.sh in.tmp $2
done > out.$1.$2.$3
fi 


# to be formatted
if [ $s == 1 ]; then
echo "# formatting set"
for j in ` awk '{print $1}' $1 |  sed 's/NM_/NM/g' |  sed 's/NP_/NP/g' | sed 's/XM_/XM/g'  | sed 's/_/ /g' | awk '{print $1}' | sort -u | sed 's/NM/NM_/g' | sed 's/NP/NP_/g'| sed 's/XM/XM_/g' `; 
do 
grep $j $1 | sort -k 2 |  uniq > in.tmp;  
sh start.sh in.tmp $2
done > out.$1.$2.$3

fi 
# computes performances at a glance (can be removed)
awk '($2>0){p++} ($2<0){n++} END{print "#", p/NR, n/NR}' out.$1.$2.$3
# reports (can be removed)
cat  out.$1.$2.$3
rm in.tmp
