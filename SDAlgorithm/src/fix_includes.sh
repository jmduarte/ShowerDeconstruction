for file in *.cc
do
  for include in ../interface/*.h    
  do
    export rule=`basename $include`
    sed -i "s:\"$rule\":\"ShowerDeconstruction/SDAlgorithm/interface/$rule\":" "$file"
  done
done