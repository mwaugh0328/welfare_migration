**What is Git? How do I use it?**

I'm still not sure! But here is my attempt at explaining what I know...

Git is a version control system. The basic idea is that it is a system which will record the changes made to a file. For example, as you modify and develop code, this is a way to retain the older version of your code, track the changes, and more! I'm still not at the more part...One final point: Git is local to your computer. This is what will distinguish it from GitHub which is a code/file/data repository that can store the code and the changes made on the world wide web.

Here are some good (much better than mine) notes by the quant econ guys.

https://github.com/jstac/quantecon_nyu_2016/blob/master/lecture2/git_intro/gitnotes.md

Step one-it is very simple to download....

* Go [here](https://git-scm.com/) and follow the instructions... I don't know much so I just followed the default options.
* Open the command prompt. And enter the following commands which will assign you a username and an email. This will be relevant when interfacing with GitHUB

```
git config --global user.name "username"    
git config --global user.email "username@email.com"
```

----------

**How to interface with GitHub**

I've been using GitHub for awhile as a place to store polished code. Below, let me talk through about how to work with git, GitHub, and then files on your computer. First, what I like to do is to create a separate file called GitHub on my computer. It is from here that anything git related will be performed on. This maybe violating the spirit of the use of git, but bear with me, I do this a lot.

Then in the command prompt I change the directory so that it is operating within, I do this so...

```
cd c:/GitHub
```

Then I'm going to clone a folder to my desktop from my GitHub account. This will be my local version of the code from which I can make changes and so forth. In this case I do the following:

```
git clone https://github.com/mwaugh0328/JIE-SW-2014
```
which creates a folder called "JIE-SW-2014" and contains all the files that are posted on GitHub. Now do the following, type in

```
cd c:/GitHub/JIW-SW-2014
git status
```

where the first line changes the directory, the second line asks what is going on with the files. It should reply back, "on branch master" and "your branch is up to date with origin/master" In other words, your files are consistent as "master" file. Now do one more thing

```
git log --pretty=oneline
```

which then is going to show the different records associated with this file. Note all the number/letter combinations, these are codes that track the different commits made.

Let me do a couple of more things. I'm going to modify one of the files. Then it is worth seeing what is going on by asking about the status.

```
git status
```
It will show a message that will say to the effect...as far as git is concerned the relevant files they have is the master, but there is changed file that has not been committed.  This is the nice feature of this. You can work locally and edit at will, but until a commit is made, the master file won't be changed. This allows you to mess stuff up at will, but always return to the master branch with one command.

Ok so how do we commit. This is a key aspect of the git, we have made a change and are going to commit it to the record. To do so you go

```
git commit -am "a test"
git status
```

which commits the modified file. And then checking the status, it should show that current changes are consistent with the master, but the branch is one step ahead of the origin/master

Now I'm going to send the change back to GitHub the code repository. To do so I do

```
git push --all
```

which then pushes these files back to where they came from. Here a window may pop up asking for your username and login for GitHub.

Here is something interesting to do, on another computer pull the folder and then check the status. What do you see? It should be something like, the origin master has the thing showing "a test".

----------
**How to Change Versions**

Ok, so this is one of the more important aspects of git. Let me walk through some commands to get a sense of what is possible.

First, checkout the log and below I'm printing the output

```    
git log --pretty=oneline
9c935cbfa518847cfdc950d5d21a125f6a87a74f (HEAD -> master, origin/master, origin/HEAD) Updated Git notes
982d88729ca4848819c51bf5b5f8f692bd6c8b49 Notes on how to use git
87cd42e4753e7bf55c83ba4733236f4084000fcd First Commit
01c0b612d80294360d80519b82e28430b3878bb2 Initial commit
```
Now the key thing to note are the number/letter combinations, these are recording the different committed versions. Not lets suppose that we want to return to a previous version. In that case you use the "checkout" command.
```
git checkout 902d88
```
And then this will revert your file to the version associated with "Notes on how to use git" commit. This is nice to see what you did in the past, etc. Then if you want to go back, you just do
```
git checkout master
```
And this will then return you to the latest commit.

Ok one more thing. Suppose you want to eliminate the previous commit and return the master to some earlier point. Note, this will permanently revert back, once you do it, all is lost(?)
```
git reset --hard 902d88
```
And then if you check the log, this will take you to
```
$ git log --pretty=oneline
982d88729ca4848819c51bf5b5f8f692bd6c8b49 (HEAD -> master, origin/master, origin/HEAD) Notes on how to use git
87cd42e4753e7bf55c83ba4733236f4084000fcd First Commit
01c0b612d80294360d80519b82e28430b3878bb2 Initial commit
```
which reverts everything back to the point where the master is not 982d88.

---

**How to Add Files**

Another important aspect of git. Ok so you save a new file or create a new folder within your repository. If you do `git status` you should see a message that says something to the following effect: on master, but there are untracked files within the repository. You want git to track those files, so what you need to do is to add them. The command is simple...

```
git add .
```
which says add files and then the `.` means add all untracked files in the repository. Then check the status again. It should say something to the effect that there you are on the master but that there are changes to be committed, and these changes shown should enumerate all the files you wanted to add. Now you just commit and you are done!

Another thing that you can do is to use wildcard notation or``*``. So for example we create a folder ``test`` and we want to add everything in just that folder we do
```
git add test\*
```
which says add everything (this is what the wildcard ``*`` does) in the test folder.

**How to Get Rid of Files**

First, this is dangerous territory, so be sure to know what you are doing before hand. So sometimes you want to reorganize your repository---this may include deleting files. What are the commands to do so. So something like

```
git rm your-file-to-delete
```
will do the trick.
