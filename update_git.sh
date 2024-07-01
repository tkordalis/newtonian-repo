#!/bin/bash

# Folder to skip
SKIP_FOLDERS=("tecplot" "1_results_dat" "sol")
for folder in "${SKIP_FOLDERS[@]}"; do
    echo "$folder"
done

# Get the current branch name
BRANCH=$(git rev-parse --abbrev-ref HEAD)

# Function to check if a folder should be skipped
should_skip() {
    for folder in "${SKIP_FOLDERS[@]}"; do
        if [[ $1 == $folder* ]]; then
            return 0
        fi
    done
    return 1
}

# Function to check if a file should be skipped based on its name
should_skip_file() {
    if [[ $1 =~ Iteration_[0-9]+\.plt$ ]]; then
        return 0
    fi
    return 1
}

# Add all files, but skip specified folder and files
find . -type f | while read file; do
    # Get the folder of the file
    folder=$(dirname "$file")

    # Check if the folder or file should be skipped
    if should_skip "$folder"; then
        echo "Skipping $file"
    else
        echo "Adding $file"
        git add "$file"
    fi
done

# Prompt the user for a commit message
echo "Enter commit message:"
read commit_message

# Commit with the provided message
git commit -m "$commit_message"

# Push to the current branch
#git push origin "$BRANCH"

