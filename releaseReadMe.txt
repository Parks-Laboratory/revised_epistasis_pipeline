# Run your jobs first

1. chmod +x release.sh
2. ./release.sh 2048 1024 600 &

# first parameter is request memory 2048 = 2G
# second parameter is reauest disk 1024 = 1G
# third one is time to repeat 600 = 10min

#  & makes it run in background

# can use ctrol + c to exit but still run in background

3. jobs  
# can see whether process is running 
4. to kill process:  kill $(jobs -p) 

