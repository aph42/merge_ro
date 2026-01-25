
for y in $(seq -f "%04g" 2020 2020)
do
   let c=$y-2007+64
   ./grid_all.sh $y 2 4 >& logs/g.e5.eqrad.${y}a.test.out &
   ./grid_all.sh $y 5 7 >& logs/g.e5.eqrad.${y}b.test.out &
   ./grid_all.sh $y 8 10 >& logs/g.e5.eqrad.${y}c.test.out &
   ./grid_all.sh $y 11 11 >& logs/g.e5.eqrad.${y}d.test.out &
   #taskset -c $c ./grid_all.sh $y 01 02 >& logs/gridtrop.${y}.T24h.out &
   #sbatch -N 1 --job-name="gridtrop.${y}.m9to10" --output="logs/gridtrop.${y}.m9to10.out" grid_all.sh $y 09 10
   #sbatch -N 1 --job-name="gridtrop.${y}.m6"  --output="logs/gridtrop.${y}.m06.out" grid_all.sh $y 06 06
   #sbatch -N 1 --job-name="gridtrop.${y}.m8"  --output="logs/gridtrop.${y}.m08.out" grid_all.sh $y 08 08
   #sbatch -N 1 --job-name="gridtrop.${y}.m10" --output="logs/gridtrop.${y}.m10.out" grid_all.sh $y 10 10
   #>& logs/grid.eq3.${y}.out &
done
