# to compile
#gcc -O3 -lm -lfann  ./scripts/one.c  -o scripts/test -L/sw/lib -lpcre -I/sw/include
# and for the path check as in
# http://stackoverflow.com/questions/4825652/how-do-i-add-a-directory-to-c-header-include-path
# pcre-config --libs --cflags

# cleans

echo "del" > new.run; rm new.run

# input file
file=$1;

cp $file new.run
sh scripts/run.sh new run > tmp/input.2.tmp
awk '(NF>2){print $0, $1, -1}' tmp/input.2.tmp > tmp/data.tmp



f="params.matrices"
# checks the length and changes the format  (here I introduced a change to control that the size "106" is correct)
awk '(NF==106)' tmp/data.tmp | awk '{for(i=1;i<=NF-2;i++){printf "%s\t",$i;} printf "%4.2f\t", 1; s=$NF; if(s<-1){s=-1} if(s>1){s=1} printf "\n%s\t%4.2f\n", $(NF-1),s}' > tmp/input.tmp
s=`wc tmp/input.tmp | awk '{print $1}'`

# if the file is non-null, then it normalizes and makes the calculations
awk '(NF>2){s=0; for(i=2;i<NF-1;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF-1;i++){if(s==0){v=a[i]} if(s!=0){v=a[i]/s} printf "%s\t", v} printf "%s\t",$NF; printf "\n"} (NF==2){print $2}'  tmp/input.tmp > tmp/input.txt

# normalizes
cat $f/min.relative.i.txt $f/max.relative.i.txt  $f/min.relative.o.txt $f/max.relative.o.txt tmp/input.txt | awk '(NR==1){for(i=1;i<=NF;i++){m[i]=$i}} (NR==2){for(i=1;i<=NF;i++){M[i]=$i}} (NR==3)&&(NF==1){m2[1]=$1} (NR==4)&&(NF==1){M2[1]=$1} (NR>4)&&(NF>1){for(i=1;i<=NF;i++){ s=0; if (m[i]!=M[i]){s=(($i-m[i])/(M[i]-m[i])-1/2)*2;} printf "%4.2f\t",s} printf "\n"} (NR>4)&&(NF==1){s=0; if (m2[1]!=M2[1]){s=(($1-m2[1])/(M2[1]-m2[1])-1/2)*2} print s}'  > tmp/norm.tmp


# writes input in the correct format

awk '{a[NR]=$0} (NF!=1){s=NF} END{printf "%i\t%i\t%i\n",NR/2,s, 1; for(i=1;i<=NR;i++){printf "%s\n",a[i]}}' tmp/norm.tmp  > input.txt


./scripts/test input.txt ./models/secondary.41 | awk '(NF==3){print $NF}'
