Data science tools
==================

UNIX/Linux operating systems: why do we care?
---------------------------------------------
- Have better facilities for automation (shell scripting)
- A lot of scientific software is (only) written for them
- Linux is free, so analytical environments can be copied
  
UNIX file system conventions
----------------------------
UNIX-like file systems:
- are hierarchical trees navigated with relative and absolute paths
- follow some organizational conventions
- consider everything a “file”
- consider files simple byte arrays (and text)

| Directory | Description          |
|-----------|----------------------|
|`/bin`     | Single user binaries |
|`/dev`     | Essential devices    |
|`/etc`     | Configuration files  |
|`/home`    | User home dirs       |
|`/lib`     | Libraries            |
|`/sbin`    | System binaries      |
|`/tmp`     | Temporary files      |
|`/usr`     | Multi user utilities |
|`/var`     | Variable files       |

File access
-----------
- Access to files is granted at user and group level
- Files can be read (r), written (w) and/or executed (x)
- To read from or write to a file, a file handle is opened
- File extensions don’t matter

Commonly-used file operations
-----------------------------

| Command | Function                |
|---------|-------------------------|
| `ls`    | List files              |
| `cd`    | Change directory        |
| `pwd`   | Print working directory |
| `rm`    | Remove files            |
| `cp`    | Copy files              |
| `mv`    | Rename/move files       |
| `mkdir` | Make directory          |
| `rmdir` | Remove directory        |
| `chmod` | Change file mode        |
| `chown` | Change file owner       |

Commonly-used text operations
-----------------------------

| Command | Function                |
|---------|-------------------------|
| `more`  | Page through text       |
| `head`  | View first lines        |
| `tail`  | View last lines         |
| `grep`  | Search for pattern      |
| `sed`   | Search and replace      |
| `sort`  | Sort lines              |
| `uniq`  | Filter duplicate lines  |
| `cat`   | Print file contents     |
| `wc`    | Word count              |
| `awk`   | Tabular data processing |

UNIX pipes
----------
- UNIX programs by default write to STDOUT and read from STDIN
- The defaults can be re-directed with “> outfile” and “< infile”, respectively
- The output from one program can be piped into the input from the next:
   
    cat file | sort | uniq > sorted_no_dups

UNIX programming
----------------
- UNIX provides a standardized environment for developing applications (historically in C)
- Many pieces of re-useable code (headers and libraries)
- Standard powerful text editors (vi, emacs)
- Standard toolchain for building binaries (e.g. autotools, make, gcc, ld, etc.)
- Standard support for scripting languages (perl, python, ruby, etc.)

Compiling
---------
Source code can be compiled into binary following standard (optional) steps:
    tar -xzf prog.tar.gz
    cd prog
    (autoconf)
    (./configure)
    make
    sudo make install

