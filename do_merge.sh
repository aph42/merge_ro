
for y in $(seq -f "%04g" 2020 2020)
do
   #nohup nice ./merge_all.sh $y 1 1 >& logs/meitest.${y}.out &
   for d in $(seq 1 41 366)
   do
      d2=$((d + 40))
      nohup nice ./merge_all.sh $y $d $d2 >& logs/m.${y}.${d}-${d2}.out &
   done
   #nohup nice ./merge_all.sh $y 1 121 >& logs/m.${y}.d1.121.out &
   #nohup nice ./merge_all.sh $y 122 244 >& logs/m.${y}.d122.244.out &
   #nohup nice ./merge_all.sh $y 245 366 >& logs/m.${y}.d245.366.out &
done
