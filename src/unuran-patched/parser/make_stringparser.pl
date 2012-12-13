#!/usr/bin/perl
##############################################################################
#                                                                            #
#           UNURAN -- Universal Non-Uniform Random number generator          #
#                                                                            #
##############################################################################
#                                                                            #
#   FILE:    make_stringparser.pl                                            #
#                                                                            #
#   Read all UNU.RAN header files and create switch tables for stringparser  #
#   and write doc file for keywords                                          #
#                                                                            #
##############################################################################
#                                                                            #
#   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold              #
#   Department of Statistics and Mathematics, WU Wien, Austria               #
#                                                                            #
#   This program is free software; you can redistribute it and/or modify     #
#   it under the terms of the GNU General Public License as published by     #
#   the Free Software Foundation; either version 2 of the License, or        #
#   (at your option) any later version.                                      #
#                                                                            #
#   This program is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#   GNU General Public License for more details.                             #
#                                                                            #
#   You should have received a copy of the GNU General Public License        #
#   along with this program; if not, write to the                            #
#   Free Software Foundation, Inc.,                                          #
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                   #
#                                                                            #
##############################################################################

use strict;

my $VERBOSE = 1;

# ----------------------------------------------------------------------------
# $Id: make_stringparser.pl 5330 2011-04-19 10:50:18Z leydold $
# ----------------------------------------------------------------------------
#
# File: make_stringparser.pl
#
##############################################################################

sub usage {
    my $progname = $0;
    $progname =~ s#^.*/##g;
        
    print STDERR <<EOM;
usage: $progname < template > Cfile

    Generates a stringparser. A template file is read from STDIN.
    It contains the auxilliary routines and the skeleton for the 
    string interpreter. The routine reads all relevant information
    from the corresponding header files and inserts the different
    cases. The output is written on STDOUT.
    Additionally the list of all key words produced for the manual
    texi format.

    Notice: The keyword _orderstatistics_ is not treated automatically
    but entered manually into the template file and in this script!
      
EOM

    exit;
}

##############################################################################
# Supported Distribution types
#
my %SUPPORTED_DISTR_TYPES =
    ( 'cont'  => 1,
      'cemp'  => 1,
      'discr' => 1,
      '-x-all-x-' => 1     # set calls unur_distr_set_... for all types
    );

# Unsupported types:
#   corder, cvec, cvemp, matr

# Commands substituted by string parser
my %SUBST_COMMANDS =
    ( 'unur_distr_cont_set_pdfstr'     => 'pdf',
      'unur_distr_cont_set_logpdfstr'  => 'logpdf',
      'unur_distr_cont_set_cdfstr'     => 'cdf',
      'unur_distr_cont_set_logcdfstr'  => 'logcdf',
      'unur_distr_cont_set_hrstr'      => 'hr',
      'unur_distr_discr_set_pmfstr'    => 'pmf',
      'unur_distr_discr_set_cdfstr'    => 'cdf' ,
      'unur_distr_discr_set_logcdfstr' => 'logcdf'
    );

# Methods ignored by string parser
my %IGNORED_METHODS =
    ( 'cext' => 1,
      'dext' => 1,
      'mixt' => 1,
      );

# Commands ignored by string parser
my %IGNORED_COMMANDS =
    ( 'unur_distr_set_extobj'       => 1,
      'unur_distr_cont_set_pdf'     => 1,
      'unur_distr_cont_set_dpdf'    => 1,
      'unur_distr_cont_set_logpdf'  => 1,
      'unur_distr_cont_set_dlogpdf' => 1,
      'unur_distr_cont_set_cdf'     => 1,
      'unur_distr_cont_set_logcdf'  => 1,
      'unur_distr_cont_set_hr'      => 1,
      'unur_distr_discr_set_pmf'    => 1,
      'unur_distr_discr_set_cdf'    => 1,
      'unur_distr_discr_set_logcdf' => 1,
      'unur_mixt_set_useinversion'  => 1,
    );

# distributions ignored by string parser
# (We cannot handle this distributions yet. Thus we will ignore them 
# until we find some time to fix this.)
my %IGNORED_DISTRIBUTIONS =
    ( 'multinormal'  => 1,   
      'multicauchy'  => 1,
      'multiexponential' => 1,
      'multistudent' => 1,
      'copula'       => 1,
      'correlation'  => 1
    );

##############################################################################

# ----------------------------------------------------------------
# Directory with sources

my $top_srcdir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';
$top_srcdir .= "/../..";

# ----------------------------------------------------------------
# Load routines for reading data about PDF from UNU.RAN files
  
require "$top_srcdir/scripts/read_PDF.pl";

# ----------------------------------------------------------------
# doc files to be generated

my $doc_dir = $ENV{'srcdir'} ? $ENV{'srcdir'} : '.';

my $doc_file = "$doc_dir/stringparser_doc.dh";

my $distr_doc_string;
my $method_doc_string;

##############################################################################
# Get all header files in methods directory
#
my $methods_dir = "$top_srcdir/src/methods";
opendir (METHDIR, "$methods_dir") or die "can't open directory $methods_dir";
my @methods_h_files = grep {/[^\#].*[.]h$/ } readdir METHDIR;
closedir METHDIR;

##############################################################################
# Get all header files in distr directory
#
my $distr_dir = "$top_srcdir/src/distr";
opendir (DISTRDIR, "$distr_dir") or die "can't open directory $distr_dir";
my @distr_h_files = grep {/[^\#].*[.]h$/ } readdir DISTRDIR;
closedir DISTRDIR;

foreach my $h (@distr_h_files) {
    print STDERR "'$h' ";
}
print STDERR "\n";

##############################################################################
# Global variables

# Store unsupported, substituted and ignored set calls for output on screen
my $msg_unsupported;
my $msg_substituted;
my $msg_ignored;
my $msg_ignored_methods;

##############################################################################
# Read template C file from STDIN and insert C code for string interpreter 
#
while ( <STDIN> ){

    unless (/^s*=INPUT\s*(\w*)/){
	print $_;
    }

    else {
	my $type = $1;    # INPUT type 

	if ( $type eq "list_of_distributions" ){
	    print make_list_of_distributions();
	}
	elsif ( $type eq "list_of_distr_sets" ){
	    print make_list_of_distr_sets();
	}
	elsif ( $type eq "list_of_methods" ){
	    print make_list_of_methods();
	}
	elsif ( $type eq "list_of_par_sets" ){
	    print make_list_of_par_sets();
	}
	else{
	    die "Error: unknown qualifier after =INPUT: $type\n";
	}
    }
}

if ($msg_substituted) {
    print STDERR "Substituted set commands:\n";
    print STDERR "$msg_substituted\n";
}

if ($msg_ignored) {
    print STDERR "Ignored set commands:\n";
    print STDERR "$msg_ignored\n";
}

if ($msg_ignored_methods) {
    print STDERR "Ignored methods:\n";
    print STDERR "$msg_ignored_methods\n";
}

if ($msg_unsupported) {
    print STDERR "Unsupported set commands:\n";
    print STDERR "$msg_unsupported\n";
}

# Print documentation.
open DOC, ">$doc_file" or die "Cannot open file $doc_file for writing";
print DOC
    "/*\n",
    "=NODE  KeysDistr   Keys for Distribution String\n",
    "=UP StringDistr [10]\n\n",
    "=DESCRIPTION\n\n",
    $distr_doc_string,
    "\n=EON\n*/\n";
print DOC
    "/*\n",
    "=NODE  KeysMethod   Keys for Method String\n",
    "=UP StringMethod [10]\n\n",
    "=DESCRIPTION\n\n",
    $method_doc_string,
    "\n=EON\n*/\n";
close DOC;

##############################################################################
#
# The end
#
exit (0);

##############################################################################

##############################################################################
#                                                                            #
# Subroutines                                                                #
#                                                                            #
##############################################################################

##############################################################################
#
# Make subroutine for getting parameter object for method
#
sub make_list_of_distributions {

    my $code;

    # List of distributions
    # For description of data fields in this list see file `read_PDF.pl'.
    my $DISTR = read_PDFdata( $top_srcdir );

    # print info on screen
    print STDERR "\nDistributions:\n" if $VERBOSE;

    # print docu
    $distr_doc_string .= "List of standard distributions "
	."\@pxref{Stddist,,Standard distributions}\n\n";
    $distr_doc_string .= "\@itemize \@minus\n";

    # make switch for first letter of distribution name
    $code .= "\t switch (*distribution) {\n";

    my $last_char;

    # Make list of all distributions
    sub caseinsensitive { ("\L$a" cmp "\L$b") };
    foreach my $distr (sort caseinsensitive keys %{$DISTR}) {

	# check whether command should be ignored
	if ($IGNORED_DISTRIBUTIONS{$distr}) {
	    # ignore this distribution
	    $msg_ignored .= "  unur_distr_$distr()\n";
	    next;
	}

	print STDERR $distr,"  ";

	my $char = lc substr $distr,0,1;

	if ($char ne $last_char) {
	    $code .= "\t\t break;\n" if $last_char;
	    $code .= "\t case '$char':\n";
	    $last_char = $char;
	}

	# print code
	$code .= "\t\t if ( !strcmp( distribution, \"\L$distr\") ) {\n";
	$code .= "\t\t\t distr = unur_distr_$distr (darray,n_darray);\n";
	$code .= "\t\t\t break;\n";
	$code .= "\t\t }\n";

	# print docu
	$distr_doc_string .= "\@item \@code{[distr =] $distr(\@dots{})} \@ \@ \@ \@ " 
	    ." \@result{} \@pxref{$distr}\n";
    }

    # end of switch for first letter
    $code .= "\t }\n\n";

    # print docu
    $distr_doc_string .= "\@end itemize\n\n\@sp 1\n";

    print STDERR "\n\n";

    # Make list of generic distribution objects
    print STDERR "Generic distributions:\n";

    # print docu
    $distr_doc_string .= "List of generic distributions "
	."\@pxref{Distribution_objects,,Handling Distribution Objects}\n\n";
    $distr_doc_string .= "\@itemize \@minus\n";

    $code .= "\t /* get pointer to generic distribution object */\n";
    $code .= "\t if (distr == (struct unur_distr *) &distr_unknown) { \n";
    $code .= "\t\t do {\n";

    foreach my $hfile (sort @distr_h_files) {
	# Read content of header file
	open H, "< $distr_dir/$hfile" or  die ("can't open file: $distr_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;
	# there must be a unur_distr_...._new call
	next unless $content =~ /[^\n]\s*unur\_distr\_(\w+)\_new/;
	# ID for method
	my $distr_type = "\L$1";
	# not all generic distributions are supported yet
	next unless $SUPPORTED_DISTR_TYPES{$distr_type};

	# make code
	print STDERR "  \U$distr_type" if $VERBOSE;
	$code .= "\t\t\t if ( !strcmp( distribution, \"$distr_type\") ) {\n";
	$code .= "\t\t\t\t distr = unur\_distr\_$distr_type\_new();\n";
	$code .= "\t\t\t\t break;\n";
	$code .= "\t\t\t }\n";

	# print docu
	$distr_doc_string .= "\@item \@code{[distr =] $distr_type} \@ \@ \@ \@ " 
	    ." \@result{} \@pxref{\U$distr_type}\n";
    }

    $code .= "\t\t } while (0);\n";
    $code .= "\t }\n\n";

    # print docu
    $distr_doc_string .= "\@end itemize\n\n\@sp 1\n";

    # add comment for order statistics ...
    $distr_doc_string .= comment_for_corder();

    # end
    print STDERR "\n\n";

    # Return result
    return $code;

} # end of make_list_of_distributions()

##############################################################################
#
# Make subroutine for setting parameters in distribution objects
#
sub make_list_of_distr_sets {

    my $set_commands;
    my $set_doc;
    my $code_unsupported;
    my $code_substituted;
    my $code_ignored;
    my $code;

    # print info on screen
    print STDERR "Set commands for Distributions:\n" if $VERBOSE;

    # Read all header files 
    foreach my $hfile (sort @distr_h_files) {

	# Read content of header file
	open H, "< $distr_dir/$hfile" or  die ("can't open file: $distr_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;

	# there must be a unur_distr_..._set call
	next unless $content =~ /[^\n]\s*unur\_distr\_.*\_set\_/;

	# distribution type
	my $distr_type;
	if ($content =~ /[^\n]\s*unur\_distr\_(\w+)\_new/) {
	    # set calls for special distribution type
	    $distr_type = "\L$1";
	}
	else {
	    # set calls for all distribution types
	    $distr_type = "-x-all-x-";
	}

	# not all generic distributions are supported yet
	next unless $SUPPORTED_DISTR_TYPES{$distr_type};
	print STDERR "  \U$distr_type: " if $VERBOSE;

	# remove obsolete functions
	$content =~ s {/\*\s*=OBSOLETE.*$} []gsx;

	# Remove all comments and empty lines ...
	$content =~ s {/\*.*?\*/} []gsx;
	$content =~ s /\n\s*\n/\n/gsx;

	# Split into lines ...
	my @lines = split /\n/, $content;

	# Get all set calls
	foreach my $l (@lines) {
	    next unless $l =~ /^\s*(\w*\s+)unur\_distr\_($distr_type\_)?set_(\w+)\s*\((.+)([^\s])\s*$/; 

	    # short name of set command
	    my $command = $3;

	    # full name of command
	    my $command_name = "unur\_distr\_".$2."set_$command";

	    # list of arguments
	    my $args = $4;

	    # Check syntax of set command
	    if ( $5 ne ';' ) {
		# last character must be a ';'
		die "Unknown syntax (terminating ';' missing) in $hfile:\n$l\n";
	    }
	    if ( $1 !~ /^int\s+$/ ) {
		# type must be 'int'
		die "Unknown syntax (function type) in $hfile:\n$l\n";
	    }
	    if ( unmatched_parenthesis($l) ) {
		# parenthesis must match
		die "Unknown syntax (umatched parenthesis) in $hfile:\n$l\n";
	    }

	    # print name of parameter
	    print STDERR "$command " if $VERBOSE;

	    # process list of args
	    $args =~ s/\)\s*$//;             # remove closing parenthesis
	    my @args_list = split /\,/, $args;

	    # first argument must be of type UNUR_DISTR
	    my $a = shift @args_list;
	    unless ($a =~ /UNUR_DISTR/) {
		die "Unknown syntax (first argument not of type UNUR_DISTR) in $hfile:\n$l\n";
	    }

	    # number of arguments
	    my $n_args = $#args_list+1;

	    # get type of arguments
	    my $type_args = "";
	    foreach my $a (@args_list) {
		my $t;
		# type of argument
		if    ($a =~ /double/)   { $t = 'd'; }
		elsif ($a =~ /int/)      { $t = 'i'; }
		elsif ($a =~ /unsigned/) { $t = 'u'; }
		elsif ($a =~ /char/)     { $t = 'c'; }
		else                     { $t = '?'; }

		# arrays are indicated by capital letters
		if ($a =~ /\*/ or $a =~ /\[/) {
		    # this is interprated as list
		    $t = "\U$t";
		}
		$type_args .= $t;
	    }

	    # check whether command should be ignored
	    if ($IGNORED_COMMANDS{$command_name}) {
		# ignore this set command
		$code_ignored .= "\t /* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		$msg_ignored .= "  $command_name()\n";
		next;
	    }

	    # make set calls
	    my $set;
	    my $doc;

	    # beginning of case
	    $set = "\t\t\t\t /* n = $n_args; type = $type_args: $args*/\n";

	    # use keyword "void" when no argument is required
	    $type_args = "void" if $type_args eq "";

	    # we support the following cases:
	    #   "i"    ... one argument of type int required
	    #   "ii"   ... two arguments of type int required
	    #   "d"    ... one argument of type double required 
	    #   "dd"   ... two arguments of type double required 
	    #   "Di"   ... a list of doubles and one argument of type int required
	    #              (the second argument is considered as size of the double array)
	    #   "C"    ... one string (array of char)

	    my %type_args_doc = 
		( 'i'  => '[= @i{<int>}]',
		  'ii' => '= @i{<int>}, @i{<int>} | (@i{<list>})',
		  'd'  => '= @i{<double>}',
		  'dd' => '= @i{<double>}, @i{<double>} | (@i{<list>})',
		  'Di' => '= (@i{<list>}) [, @i{<int>}]',
		  'C'  => '= "@i{<string>}"'
		  );

	    if ($type_args =~ /^(i|ii|d|dd|Di|C)$/) {
		my $type = $1;
		$set .= "\t\t\t\t result = _unur_str_distr_set_$type(distr,key,type_args,args,$command_name);\n";
		unless ($set_commands->{$distr_type}->{$command}) {
		    $set_commands->{$distr_type}->{$command} = $set; }
		else {
		    die "\nset command redefined: $distr_type/$command"; }

		# check whether command should also have substitute
		if ($SUBST_COMMANDS{$command_name}) {
		    my $command_subst = $SUBST_COMMANDS{$command_name};
		    $msg_substituted .= "  $command_name()  --> $SUBST_COMMANDS{$command_name}\n";
		    $code_substituted .= "\t /* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		    unless ($set_commands->{$distr_type}->{$command_subst}) {
			$set_commands->{$distr_type}->{$command_subst} = $set; }
		    else {
			die "\nset command redefined: $distr_type/$command_subst"; }
		}

		# make docu
		if ($SUBST_COMMANDS{$command_name}) {
		    $command = $SUBST_COMMANDS{$command_name}; }
		$set_doc->{$distr_type}->{$command} = 
		    "\@item $command $type_args_doc{$type_args}\n \@result{} "
		    ."\@pxref{funct:$command_name,,\@command{$command_name}}\n";
	    }

	    else {
		# cannot handle this set command
		$code_unsupported .= "\t/* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		$msg_unsupported .= "  $command_name()\n";
	    }
	}

	# end of distribution type
	print STDERR "\n" if $VERBOSE;
    }

    # print info on screen
    print STDERR "\n" if $VERBOSE;


    # get list of all distribution types
    my @distr_type_list = sort (keys %{$set_commands});

    # print docu
    $distr_doc_string .= "List of keys that are available via the String API.\n"
	."For description see the corresponding UNU.RAN set calls.\n\n";
    $distr_doc_string .= "\@itemize \@bullet\n";

    # first distribution type MUST be -x-all-x-
    my $dt = shift @distr_type_list;
    unless ($dt eq "-x-all-x-") {
	die "first distribution type MUST be '-x-all-x-'\n";
    }
    # List of set calls for all distribution types
    $distr_doc_string .= "\@item All distribution types\n";
    my $code_x_all_x_ = make_distr_set_calls($dt,$set_commands,$set_doc);

    # List of set calls for distribution types
    # switch for distribution types
    $code .= "\n\t switch (distr->type) {\n";
    foreach $dt (@distr_type_list) {

	# print docu
	$distr_doc_string .= "\@item \@code{$dt} \@ \@i{(Distribution Type)}\@ \@ \@ \@ "
		."(\@pxref{\U$dt})\n";

	# make label for distribution type
	$code .= "\t case UNUR_DISTR_\U$dt:\n";

	# make list of all set commands for distribution type
	$code .= make_distr_set_calls($dt,$set_commands,$set_doc);

	# end of case for distribution type
	$code .= "\t\t break;\n";
    }
    # end of switch for distribution types
    $code .= "\t }\n";

    # append list for set calls for all distribution types
    $code .= "\n\t /* set calls for all distribution types */\n";
    $code .= "\t if (result == UNUR_ERR_STR_UNKNOWN) {\n";
    $code .= $code_x_all_x_;
    $code .= "\t }\n\n";

    # print docu
    $distr_doc_string .= "\@end itemize\n\n";

    # add comment on igored and unsupported code into C file
    if ($code_ignored) {
	$code .= "\n\t /* Ignored set commands: */\n $code_ignored\n"; }
    if ($code_substituted) {
	$code .= "\n\t /* Subsituted set commands: */\n $code_substituted\n"; }
    if ($code_unsupported) {
	$code .= "\n\t /* Unsupported set commands: */\n $code_unsupported\n"; }

    # Return result
    return $code;

} # end of make_list_of_distr_sets() 

##############################################################################
#
# Make subroutine for getting parameter object for method.
# Print set command.
#
sub make_distr_set_calls { 
    my $dt = $_[0];
    my $set_commands = $_[1];
    my $set_doc = $_[2];

    my $code;

    # make list of set commands for distributions
    my @command_list = sort (keys %{$set_commands->{$dt}});
    
    # make switch for first letter of key name
    $code .= "\t\t switch (*key) {\n";

    # print docu
    $distr_doc_string .= "\@table \@code\n";
    
    my $last_char;
    
    foreach my $c (@command_list) {
	
	my $char = substr $c,0,1;
	
	if ($char ne $last_char) {
	    $code .= "\t\t\t break;\n" if $last_char;
	    $code .= "\t\t case '$char':\n";
	    $last_char = $char;
	}
	
	$code .= "\t\t\t if ( !strcmp(key, \"$c\") ) {\n";
	$code .= $set_commands->{$dt}->{$c};
	$code .= "\t\t\t\t break;\n";
	$code .= "\t\t\t }\n";
	
	# print docu
	$distr_doc_string .= $set_doc->{$dt}->{$c};
    }
    
    # end of switch for first letter
    $code .= "\t\t }\n";

    
    # add comment about order statistics
    if ($dt eq "cont") {
	$distr_doc_string .= comment_for_orderstatistics();
    }

    # end of table
    $distr_doc_string .= "\@end table\n\n\@sp 1\n";

    # Return result
    return $code;

} # end of make_distr_set_calls()

##############################################################################
#
# Make subroutine for getting parameter object for method
#
sub make_list_of_methods {

    my $code;

    # Get list of all methods
    my @method_list;

    # Read all header files
    foreach my $hfile (sort @methods_h_files) {

	# Read content of header file
	open H, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;

	# We skip over all header files that do not correspond
	# to a method.
	next unless $content =~ /[^\n]\s*=METHOD\s+(\w+)/;

	# save ID for method
	push @method_list, "\L$1";
    }

    # sort list of methods
    @method_list = sort @method_list;

    # make code 

    # make switch for first letter of method name
    $code .= "\t switch (*method) {\n";

    my $last_char;

    foreach my $method (@method_list) {

	# check whether method should be ignored
	if ($IGNORED_METHODS{$method}) {
	    # ignore this method
	    $msg_ignored_methods .= "  unur_$method\_new()\n";
	    next;
	}

	my $char = substr $method,0,1;

	if ($char ne $last_char) {
	    $code .= "\t\t break;\n" if $last_char;
	    $code .= "\t case '$char':\n";
	    $last_char = $char;
	}

	# print code
	$code .= "\t\t if ( !strcmp( method, \"$method\") ) {\n";
	$code .= "\t\t\t par = unur_$method\_new(distr);\n";
	$code .= "\t\t\t break;\n";
	$code .= "\t\t }\n";
    }

    # end of switch for first letter
    $code .= "\t }\n";

    # Return result
    return $code;

} # end of make_list_of_methods()

##############################################################################
#
# Make subroutine for setting parameters in parameter objects
#
sub make_list_of_par_sets {

    my $set_commands;
    my $set_doc;
    my $code_unsupported;
    my $code_substituted;
    my $code_ignored;
    my $code;

    # print info on screen
    print STDERR "Set commands for Methods:\n" if $VERBOSE;

    # Read all header files 
    foreach my $hfile (sort @methods_h_files) {

	# Read content of header file
	open H, "< $methods_dir/$hfile" or  die ("can't open file: $methods_dir/$hfile");
	my $content = '';
	while (<H>) { $content .= $_; } 
	close H;

	# We skip over all header files that do not correspond
	# to a method.
	next unless $content =~ /[^\n]\s*=METHOD\s+(\w+)/;

	# ID for method
	print STDERR "  $1: " if $VERBOSE;
	my $method = "\L$1";

	# remove obsolete functions
	$content =~ s {/\*\s*=OBSOLETE.*$} []gsx;

	# Remove all comments and empty lines ...
	$content =~ s {/\*.*?\*/} []gsx;
	$content =~ s /\n\s*\n/\n/gsx;

	# Split into lines ...
	my @lines = split /\n/, $content;

	# Get all set calls
	foreach my $l (@lines) {
	    next unless $l =~ /^\s*(\w*\s+)unur_$method\_set_(\w+)\s*\((.+)([^\s])\s*$/; 

	    # short name of set command
	    my $command = $2;

	    # full name of command
	    my $command_name = "unur\_$method\_set_$command";

	    # list of arguments
	    my $args = $3;

	    # Check syntax of set command
	    if ( $4 ne ';' ) {
		# last character must be a ';'
		die "Unknown syntax (terminating ';' missing) in $hfile:\n$l\n";
	    }
	    if ( $1 !~ /^int\s+$/ ) {
		# type must be 'int'
		die "Unknown syntax (function type) in $hfile:\n$l\n";
	    }
	    if ( unmatched_parenthesis($l) ) {
		# parenthesis must match
		die "Unknown syntax (umatched parenthesis) in $hfile:\n$l\n";
	    }

	    # print name of parameter
	    print STDERR "$command " if $VERBOSE;

	    # process list of args
	    $args =~ s/\)\s*$//;             # remove closing parenthesis
	    my @args_list = split /\,/, $args;

	    # first argument must be of type UNUR_PAR
	    my $a = shift @args_list;
	    unless ($a =~ /UNUR_PAR/) {
		die "Unknown syntax (first argument not of type UNUR_PAR) in $hfile:\n$l\n";
	    }

	    # number of arguments
	    my $n_args = $#args_list+1;

	    # get type of arguments
	    my $type_args = "";
	    foreach my $a (@args_list) {
		my $t;
		# type of argument
		if    ($a =~ /double/)   { $t = 'd'; }
		elsif ($a =~ /int/)      { $t = 'i'; }
		elsif ($a =~ /unsigned/) { $t = 'u'; }
		elsif ($a =~ /char/)     { $t = 'c'; }
		else                     { $t = '?'; }

		# arrays are indicated by capital letters
		if ($a =~ /\*/ or $a =~ /\[/) {
		    # this is interprated as list
		    $t = "\U$t";
		}
		$type_args .= $t;
	    }

	    # check whether command should be ignored
	    if ($IGNORED_COMMANDS{$command_name}) {
		# ignore this set command
		$code_ignored .= "\t /* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		$msg_ignored .= "  $command_name()\n";
		next;
	    }
	    # check whether method should be ignored
	    for my $m (keys %IGNORED_METHODS) {
		if ($command_name =~ /^unur_$m/) {
		    # ignore this method
		    $msg_ignored .= "  $command_name()\n";
		    next;
		}
	    }

	    # make set calls
	    my $set;

	    # beginning of case
	    $set = "\t\t\t\t /* n = $n_args; type = $type_args: $args*/\n";

	    # use keyword "void" when no argument is required
	    $type_args = "void" if $type_args eq "";

	    # we support the following cases:
	    #   void   ... no argument required
	    #   "i"    ... one argument of type int required
	    #   "ii"   ... two arguments of type int required
	    #   "u"    ... one argument of type unsigned required
	    #   "d"    ... one argument of type double required 
	    #   "dd"   ... two arguments of type double required 
	    #   "iD"   ... one argument of type int and a list of doubles required
	    #              (the first argument is considered as size of the double array)
	    #   "Di"   ... a list of doubles and one argument of type int required

	    my %type_args_doc = 
		( 'void' => ' ',
		  'i'  => '[= @i{<int>}]',
		  'ii' => '= @i{<int>}, @i{<int>} | (@i{<list>})',
		  'u'  => '= @i{<unsigned>}',
		  'd'  => '= @i{<double>}',
		  'dd' => '= @i{<double>}, @i{<double>} | (@i{<list>})',
		  'iD' => '= @i{<int>} [, (@i{<list>})] | (@i{<list>})',
		  'Di' => '= (@i{<list>}), @i{<int>}'
		  );

	    if ($type_args =~ /^(void|i|ii|u|d|dd|iD|Di)$/) {
		my $type = $1;
		if ($type_args =~ /^(iD|Di)$/) {
		    $set .= "\t\t\t\t result = _unur_str_par_set_$type(par,key,type_args,args,$command_name,mlist);\n"; }
		else {
		    $set .= "\t\t\t\t result = _unur_str_par_set_$type(par,key,type_args,args,$command_name);\n"; }
		unless ($set_commands->{$method}->{$command}) {
		    $set_commands->{$method}->{$command} = $set; }
		else {
		    die "\nset command redefined: $method/$command"; }

		# check whether command should also have substitute
		if ($SUBST_COMMANDS{$command_name}) {
		    my $command_subst = $SUBST_COMMANDS{$command_name};
		    $msg_substituted .= "  $command_name()  --> $SUBST_COMMANDS{$command_name}\n";
		    $code_substituted .= "\t /* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		    unless ($set_commands->{$method}->{$command_subst}) {
			$set_commands->{$method}->{$command_subst} = $set; }
		    else {
			die "\nset command redefined: $method/$command_subst"; }
		}

		# make docu
		if ($SUBST_COMMANDS{$command_name}) {
		    $command = $SUBST_COMMANDS{$command_name}; }
		$set_doc->{$method}->{$command} = 
		    "\@item $command $type_args_doc{$type_args}\n \@result{} "
		    ."\@pxref{funct:$command_name,,\@command{$command_name}}\n";
	    }

	    else {
		# cannot handle this set command
		$code_unsupported .= "\t /* $l\n\t\t n = $n_args; type = $type_args\t */\n";
		$msg_unsupported .= "  $command_name()\n";
	    }
	}

	# end of method
	print STDERR "\n" if $VERBOSE;
    }

    # print info on screen
    print STDERR "\n" if $VERBOSE;


    # get list of all methods 
    my @method_list = sort (keys %{$set_commands});

    # make switch for methods
    $code .= "\t switch (par->method) {\n";

    # print docu
    $method_doc_string .= "List of methods and keys that are available via the String API.\n"
	."For description see the corresponding UNU.RAN set calls.\n\n";
    $method_doc_string .= "\@itemize \@bullet\n";

    foreach my $m (@method_list) {

	# print docu
	$method_doc_string .= "\@item \@code{method = $m} \@ \@ \@ \@ "
	    ." \@result{} \@command{unur\_$m\_new}\n"
	    ."(\@pxref{\U$m})\n";

	# make list of set commands for method 
	my @command_list = sort (keys %{$set_commands->{$m}});

	# make label for method
	$code .= "\t case UNUR_METH_\U$m:\n";
	
	# make switch for first letter of key name
	$code .= "\t\t switch (*key) {\n";

	# print docu
	$method_doc_string .= "\@table \@code\n";

	my $last_char;
	foreach my $c (@command_list) {

	    my $char = substr $c,0,1;

	    if ($char ne $last_char) {
		$code .= "\t\t\t break;\n" if $last_char;
		$code .= "\t\t case '$char':\n";
		$last_char = $char;
	    }

	    $code .= "\t\t\t if ( !strcmp(key, \"$c\") ) {\n";
	    $code .= $set_commands->{$m}->{$c};
	    $code .= "\t\t\t\t break;\n";
	    $code .= "\t\t\t }\n";

	    # print docu
	    $method_doc_string .= $set_doc->{$m}->{$c};
#	    $method_doc_string .= "\@item $c\n \@result{} "
#		."\@pxref{funct:unur\_$m\_set\_$c,,\@command{unur\_$m\_set\_$c}}\n";
	}

	# end of switch for first letter
	$code .= "\t\t }\n";
	$code .= "\t\t break;\n";

	# print docu
	$method_doc_string .= "\@end table\n\n\@sp 1\n";
    }

    # end of switch for methods
    $code .= "\t }\n";

    # print docu
    $method_doc_string .= "\@end itemize\n\n";

    # add comment on igored and unsupported code into C file
    if ($code_ignored) {
	$code .= "\n\t /* Ignored set commands: */\n $code_ignored\n"; }
    if ($code_substituted) {
	$code .= "\n\t /* Subsituted set commands: */\n $code_substituted\n"; }
    if ($code_unsupported) {
	$code .= "\n\t /* Unsupported set commands: */\n $code_unsupported\n"; }

    # Return result
    return $code;

} # end of make_list_of_par_sets() 

##############################################################################
#
# Simple test for unmatched parenthesis.
# It returns the number of opening parenthesis `(' minus the number
# of closing parenthesis `)' in argument $_.
#
sub unmatched_parenthesis {
    my $open = 0;   # number of unmatched `('

    ++$open while /\(/g;
    --$open while /\)/g;

    return $open;
} # end of unmachted_parenthesis()

##############################################################################
#
# Some auxiliary comments.
#

sub comment_for_corder {

    return <<EOS;
\@emph{Notice}:
Order statistics for continuous distributions (\@pxref{CORDER}) are
supported by using the key \@code{orderstatistics} for distributions
of type \@code{CONT}.

\@sp 1

EOS

} # end if comment_for_corder()

##############################################################################

sub comment_for_orderstatistics {

    return <<EOS;
\@item orderstatistics = \@i{<int>}, \@i{<int>} | (\@i{<list>})
    Make order statistics for given distribution. The first parameter
    gives the sample size, the second parameter its rank.
    (see \@pxref{funct:unur_distr_corder_new,,\@command{unur_distr_corder_new}})
EOS

} # end if comment_for_orderstatistics()

##############################################################################
