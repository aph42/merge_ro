#for y in $(seq -f "%04g" 2016 1 2018)
for y in $(seq -f "%04g" 2019 2025)
do
   #nohup nice ./proc_mission.sh $y $y >& logs/pr.geo.spire.planetiq.${y}.out &
   #nohup nice ./proc_mission.sh $y $y >& logs/pr.cosmic2021.${y}.out &
   nohup nice ./proc_mission.sh $y $y >& logs/pr.cosmic2.${y}.reproc.out &
   #nohup nice ./dl_mission.sh $y $((y + 2)) >& logs/dl.cosmic2021.${y}.out &
   #nohup nice ./dl_mission.sh $y $((y + 2)) >& logs/dl.tdx.${y}.out &
done
