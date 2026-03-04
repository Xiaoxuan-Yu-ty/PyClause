# Collection of git commands for various purposes

## Project Reconstructing
Scenario: Moving an existing git repository (eg. Pyclause) to a parent folder (eg. KGRules) as a subfolder, and add new subfolders under the parent folder as a single, unified Git repository.

1. Create physical structures
    + create a root folder: `mkdir KGRules`
    + move the old project folder into it: `mv Pyclause KGRules\`
    + create or move other folders under new root folder
2. Promote the Git logic
   + Git tracking usually lives inside your project folder in a hidden .git directory. To track everything, it must move to the top.
   + run inside new root directory (KGRules): `mv Pyclause/.git .`
3. Handling large files
   git push can fail due to commited large files, now need to recall the old commit and add large files to .gitignore and then commit
   +  Recall (reset) the commit: 
      +  `git reset --soft HEAD~1` -- reset the last commit;
      +  `git reset --mixed origin/branch_name` -- reset to the current branch's state, which can recall all failed pushes.
   +  Force un-track: `git rm -r -f --cached .`
   +  Add big files (> 100mb) to .gitignore by copy paste;
   +  Re-Add small filed to git: `git add .`
   +  Verify added changes: `git status` -- now need to look at all the changes to make sure no large files are there;
   +  Push: `git push origin branch_name`