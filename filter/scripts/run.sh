
protein=$1; val=$2; 

v=`echo $val | awk '($0~/pos/){print +1} ($0~/neg/){print -1}'`



A=-50; B=50; th=1;
#for rna in `cat $protein.$val | sort -nk 3 | sed 's/_/ /g' | awk '{print $3"_"$4}' | sort -u`;
for rna in `cat $protein.$val | sort -nk 3 | perl -ne '$_=~ m/^\S+\_\d+\-\d+\s+(\S+)_\d+\-\d+\s+\S+\s+\S+$/g;$rna=$1;print $rna,"\n";' | sort -u`;

do 
grep $rna $protein.$val  > tmp/test.txt

awk 'BEGIN{ b=int('$A'); e=int('$B');p=int('$th'); a[b-1]=0; for(i=b;i<=e;i=i+p){a[i]=0; k++} a[e+1]=0; printf "%s_%s\t", "'$protein'", "'$rna'"} 

{for(i=b;i<=e;i=i+p){min=i; max=min+p; if(($3>=min)&&($3<max)){a[i]++}} if($3<b){a[b-1]++}   if($3>e){a[e+1]++}   } 

END{printf "%i\t", a[b-1]; for(i=b;i<=e;i=i+p){printf "%i\t",a[i]}} END{printf "%i\n", a[e+1]}' tmp/test.txt; 

echo  $protein"_"$rna $v

done
