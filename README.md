# lcls-ii-geo-calculator
Python tools for diffractometer geometry calculations for LCLS-II 
Contact sarah.c.spaugh@gmail.com for implementation questions.

Functionalities include (most all still under development to some extent):
- data input via CIF file (adding vasp file support)
- calculation of goniometer/detector angles based on specified incidence or exit angle
- visualization of samples in 3d space relative to lab space

# Requirements
See requirements.txt for necessary packages.
This project is written for python 3 and in particular the GUI may not run in python 2.

# Running as a shell application
From terminal/bash shell in main directory, run main.py to launch the calculator in the terminal. 
Commands can be listed via the '?' command. Typing 'help [command name]' will print additional information for specific commands.
General syntax is [command_name] [inputs]

# Running the GUI
The file GeometryCalcsGUI will start the interface when run.

# Notes for future development
 tbd 