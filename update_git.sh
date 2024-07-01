#!/bin/bash

# Get the current branch name
BRANCH=$(git rev-parse --abbrev-ref HEAD)

git add .

# Prompt the user for a commit message
echo "Enter commit message:"
read commit_message

# Commit with the provided message
git commit -m "$commit_message"

# Push to the current branch
#git push origin "$BRANCH"

