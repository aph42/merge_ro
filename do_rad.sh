
for y in 2010
#for y in $(seq -f "%04g" 2006 2006)
do
   for d in $(seq 1 41 366)
   do
      d2=$((d + 40))
      #nohup nice ./merge_all.sh $y $d $d2 >& logs/m.${y}.${d}-${d2}.out &
      nohup nice ./rrtm_all.sh  $y $d $d2 >& logs/rad.pyr.${y}.${d}-${d2}.out &
   done
done
