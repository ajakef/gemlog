# Git Reference
All git commands must be run after navigating to the repository folder. The only exception is `git clone`, which creates the repository folder.

# Overview of gemlog editing workflow
Unless you're a registered collaborator on gemlog, you need to log in to github, fork the gemlog repository, push your work to your fork, and then send a pull request to the primary github repository.

If you are a registered collaborator on gemlog, you don't need to fork it, and you can push straight to the primary github repo.

Either way, your code should be on a named branch--please don't push code to 'main'.

# Absolute Basic Route Tasks
### Update a repository you already have
```
git pull
```
### Commit your changes to the current branch (creates a "save point")
```
git commit -a -m 'write a description of the commit changes here in quotes'
```

### Push your changes on the 'bug_fix' branch to your own gemlog fork
```
git push origin bug_fix
```

# Next Level Basic Routine Tasks

### See what branches already exist
```
git branch -v
```

### Switch to an existing branch (in this case, named 'bug_fix')
```
git checkout bug_fix
```

### Create a new branch (in this case, named 'bug_fix')
```
git checkout -b bug_fix
```

### Move your uncommitted work to existing branch bug_fix
```
git stash
git checkout bug_fix
git stash apply
```

### Start tracking a new file 'new_file.py' (must be done before committing)
```
git add new_file.py
```

# Basic One-Time Tasks

### Download the repository for the first time:
```
git clone https://github.com/ajakef/gemlog/
```

### See what remotes are already set up
```
git remote -v
```

### Set the new "upstream" remote:
This is where you download updates from. 

```
git remote add upstream https://github.com/ajakef/gemlog
```

### Change the existing "origin" remote:
This is where you upload your work to. If you're registered as a collaborator, you can push code directly to the primary repository (ajakef/gemlog). Otherwise, you can push code to your own gemlog fork on your own github site; in that case, just replace 'ajakef' below with your own github username.

```
git remote set_url origin https://github.com/ajakef/gemlog
```

# Occasional Intermediate Tasks
### Check how current the repository is
```
git status
```

### See all uncommitted changes
```
git diff
```

### Reverting all uncommitted changes
Careful with this--there's no undoing it. This command must be run from the root folder of the repository.

```
git checkout .
```

### List a log of commits leading to this point (partial sample output below)

```
git log 
commit d875a48995b2f532b85b697463aa69dc12493d95 (HEAD -> main, tag: 1.3.6, origin/main, origin/HEAD)
Author: Jake Anderson <ajakef@gmail.com>
Date:   Tue Jan 12 15:56:04 2021 -0700

    minor fix for inventory compliance with IRIS requirements

commit 8890ddeaa6b6e15af67ffacfe14e63da4fd56824
Author: Jake Anderson <ajakef@gmail.com>
Date:   Tue Jan 12 15:46:03 2021 -0700

    minor fix to gem_network

commit cee5b2e1dd5890cc2240f4e595e0415be0f3ad4a
Author: Jake Anderson <ajakef@gmail.com>
Date:   Mon Jan 11 14:51:01 2021 -0700

    minor changes to README.md and Installation.md
```

### Temporarily check out an old commit (in this case, labeled cee5b2e1dd5890cc2240f4e595e0415be0f3ad4a)
It will tell you you're in "detached HEAD" state. This means that you can see the code in the commit you're checking out, and you can modify the code create a new commit based on it, but the new commit won't belong to any branch! You probably don't want that, so if you want to start creating new code based on an old commit, you should create a new branch here before committing anything.

```
git checkout cee5b2e1dd5890cc2240f4e595e0415be0f3ad4a
```

### Return to what you were working on before checking out an old commit
This assumes you were working on branch 'bug_fix'.

```
git checkout bug_fix
```

### Amend the last commit
You can make a quick fix to the last commit if you found a small bug or if the commit message was wrong.

```
git commit -a --amend -m 'new commit message'
```


### Un-commit the last commit
This won't delete the work you did, but it will remove the commit from the log. You will need to re-commit any changes you make.

```
git reset HEAD^
```

