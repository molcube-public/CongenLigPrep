# Loop over all directories in the current directory
for dir in */; do
  # Check if the item is a directory
  if [ -d "$dir" ]; then
    # Enter the directory
    cd "$dir" || exit

    # Run your desired command here
    echo "Running something in directory: $dir"
    
    bash run-ketcher2sdf*.sh 

    # Move back to the parent directory
    cd ..

  fi
done


