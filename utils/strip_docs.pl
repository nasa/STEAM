#! /usr/bin/perl

my $root = $ARGV[0];
my $debug = 0;
   $debug = 1 if (@ARGV > 1);

print "Working on $root.C\n";

open INPUT, "$root.C";

my $func   = "";
my $docstr = "";

my %docstrs;

while (<INPUT>) {

   chomp();
   # Is this the start of a doc string?
   if (/\"\"\"/ && $docstr eq "") {
      $docstr .= "$_\n";
      next;
   }

   # Is this the end of a doc string?
   if (/\"\"\"/ && $docstr ne "") {
      $docstr .= "$_\n";
      # Save and clear
      print " DocStr for $func :\n$docstr\nDONE\n" if ($debug > 0);
      $docstrs{$func} = $docstr;
      $docstr = "";
      next;
   }
   
   # If I'm in a docstring, then just eat stuff
   if ($docstr ne "") {
      $docstr .= "$_\n";
      next;
   }

   # Look for a function name:
   my ($search) = $_ =~ /\s([:a-zA-Z0-9_-]+)\s*\(/;
   $func = $search if ($search ne "");
}
close INPUT;

# OK, so I should have them all in a hash.  Write them to the
# appropriate doc file with SWIG declarations.
open OUT, ">$root.docs";

print OUT "%feature(\"autodoc\", \"0\");\n";

my $i = 0;
foreach my $func (keys %docstrs) {
   print " Writing docstr for $func\n";
   
   # Clean it out a bit:
   $docstr = $docstrs{$func};
   $docstr =~ s/\"\"\"//g;
   $docstr =~ s/^\s*\Q*\E//g;  # Fist occurance in string
   $docstr =~ s/\n\s*\Q*\E/\n/g; # At start of each new line

   # Indentation is a bit annoying.  Keywords should not be
   # indented, but what come afterwards should be
   foreach $key (qw/Args Returns Note Raises/) {
      # Get rid of leading white spaces, but not blank lines.
      $docstr =~ s/\n[^\n]*$key.*:/\n $key:/g;
   }
   print " DocStr for $func :\n$docstr\nDONE\n" if ($debug > 0);

   $i +=1;
   my $label = "DOCSTRING$i";
   print OUT "%define $label\n\"$docstr\"\n%enddef\n";
   print OUT "%feature(\"docstring\", $label) $func;\n";
}
