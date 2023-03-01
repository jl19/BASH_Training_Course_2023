# Introduction to Linux in HPC


Linux is the best-known and most-used open source operating system.

### Historical Background

1969: Unix operating system is developed at Bell Laboratories by Ken Thompson, Dennis Ritchie and others.

* Written in C
* Already a successor to Multics
* Over time, many variation developed


1990: Linux is created  by Linus Torvalds and the Free Software Foundation.

* Interface that all Unix systems implement
* Adopted by Unix-like systems
* Perform variant calling according to GATK best practices
* Perform a variant annotation


### Popular Linux Distributions

* Debian
* Fedora and Red Hat
* Ubuntu
* Linux Mint

### Computers with Linux:

* 500 out of the top 500 supercomputers (2020)
* Web servers (95%)
* Monbile devices: 60-80% of mobile devices
* Desktop PCs: 1-2%


### The Command Line

- A line where you type commands
- Other terms:  CLI (command line interface), Console/Terminal,shell
 
Advantages

  - Simple
  - Easy to program
  - Fast and efficient to use once you know it
  
Disadvantage

  - lots of commands to remember
  
***********************************
### Navigating the Command line
***********************************

- Recalling commands

Up, down, left and right arrows of the keyboard

![keyboard1](./assets/Keyboard1.png)

![The reference diagram](./Docs/assets/navigate-command-line-1.png)
***********************************
- Completing commands/filenames

Using TAB button to auto-complete commands and filenames

![keyboard2](/Docs/assets/Keyboard2.png)

### Exercises in Terminal

- What is your username on a linux Computers

In the command line prompt, e.g.
   
```  
[username@hostname ~]$ 

[jl19@spectre10 ~]$
 
```
you can find your username as username directly. The Linux command whoami can also show your username.

```  
[username@hostname ~]$ whoami
 
```

- What's the hostname of a Linux computer
In the command line prompt, e.g.
```  
[username@hostname ~]$ hostname
```  
you can find the hostname as hostname directly. The Linux command hostname can also give you the hostname.

- The Linux command pwd can also show your current working directory
     The up- and down- (↑ and ↓) arrow keys can be used to navigate command history.

``` 
[jl19@spectre10 ~]$ pwd
/home/j/jl19
[jl19@spectre10 ~]$ 

``` 

### Directory Structure

Path: location inside file systems

``` 
[jl19@spectre10 ~]$ ls 
[jl19@spectre10 ~]$ ls /home
a  c  e  g  i  k  m  o  q  s  u  w  y
b  d  f  h  j  l  n  p  r  t  v  x  z
[jl19@spectre10 ~]$ clear

[jl19@spectre10 ~]$ ls -l/home
```

Absoute path starts with /

``` 
[jl19@spectre10 ~]$ cd
```

Relative path: relative to working directory
 
 
###Linux and NGS (Next Genearation Sequenceing) Analysis
 
Many tools developed for NGS data processing and analysis are based on Unix operating systems such as Linux. Command line tools such as awk, grep, sed are just examples of command line tools are used.

In Linux, Shell is a program that takes your commands from the keyboard and gives them to the operating system. Most Linux systems utilize Bourne Again SHell (bash), but there are several additional shell programs on a typical Linux system such as ksh, tcsh, and zsh. To see which shell you are using, type

```
[jl19@spectre10 ~]$ echo $SHELL
/bin/bash
```



### [Linux Commands Sheet](linux-command-sheets.md)

