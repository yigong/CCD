#!/usr/bin/perl
#
# This script executes a command from a list of files based on the current MPI id.
#
# Last modified: Mar/11/2005
#

# call getid to get the MPI id number

($myid,$numprocs) = split(/\s+/,`./getid`);
$file_id = $myid;
$file_id++;


# open file and execute appropriate command

$file_to_use = $ARGV[0];
open (INPUT_FILE, $file_to_use) or &showhelp;

for ($i = 1; $i <= $file_id; $i++)
{
        $buf = <INPUT_FILE>;
}

system("$buf");

close INPUT_FILE;


sub showhelp
{
        print "\nUsage: my_script.pl <filename>\n\n";
        print "<filename> should contain a list of executables, 
one-per-line, including the path.\n\n";
}
