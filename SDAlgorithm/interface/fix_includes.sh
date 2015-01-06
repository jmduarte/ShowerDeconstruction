for file in *.h
do
  for rule in *.h    
  do
    sed -i "s:\"$rule\":\"ShowerDeconstruction/SDAlgorithm/interface/$rule\":" "$file"
  done
done