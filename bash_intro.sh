#!/usr/bin/env bash

echo "Hello $(whoami)"
echo "This computer is called $(hostname)"
echo "The time is $(date)"

FOO="bar"
echo "The value of FOO is $FOO"

echo "A rainbow:"
for COL in 'red' 'orange' 'yellow' 'green' 'cheesecake' 'blue' 'indigo' 'violet'
do 
  if test $COL == 'cheesecake'; then
    echo "Cheesecake is NOT a colour of the rainbow."
  elif test $COL == 'green' ; then
    echo "$COL Bay Woooooo!"
  else
    echo $COL
  fi
done

for NUM in $(seq 0 10 100)
do
  echo $NUM
done

for FILE in /*
do
  echo $FILE
done


#interactive input
#read -p 'Enter 2 integers' ONE TWO
#echo "You entered $ONE $TWO"

#positional params
echo "First parameter is $1"
echo "All parameters are  $@"

