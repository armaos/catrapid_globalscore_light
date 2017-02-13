th=$2; 
 A=$3; 
 B=$4;
 f=$1; 
 g=`echo $f    | sed 's/cat/arr/g'`;
cm=`head -1 $f | awk '{print NF+1}'`; 

# counts the values
for c in $(seq 2 $cm); do awk 'BEGIN{a['$A'-'$th']=0; for(i='$A';i<='$B';i=i+'$th'){a[i]=0 } a['$B'+'$th']=0} (NR==1){printf "%s\t",$('$c'-1)} (NR>1){ for(i='$A';i<='$B';i=i+'$th'){min=i; max=min+'$th'; if(($'$c'>=min)&&($'$c'<max)){a[i]++}} if($'$c'<'$A'){a['$A'-'$th']++}   if($'$c'>'$B'){a['$B'+'$th']++}   } END{printf "%i\t", a['$A'-'$th']; for(i='$A';i<='$B';i=i+'$th'){printf "%i\t",a[i]}} END{printf "%i\n", a['$B'+'$th']}' $f;  done > tmp/test.$f.bis2


# assembles the file
for j in `cat tmp/test.$f.bis2 | awk '{print $1}'`; 
do 
grep $j tmp/test.$f.bis2 | head -1
grep $j $g | head -1 | awk '{print $1, (log($2)+log($3))/2}' 
done > tmp/data.$f.bis2

