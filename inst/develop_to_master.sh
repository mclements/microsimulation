#!/bin/bash
# http://stackoverflow.com/questions/2763006/change-the-current-branch-to-master-in-git

git checkout develop
git merge --strategy=ours master    # keep the content of this branch, but record a merge
git checkout master
git merge develop             # fast-forward master up to the merge
