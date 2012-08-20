#!/usr/bin/perl

# Copyright (C) 2008 by CPqD

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

BEGIN { $ENV{PATH} = '/bin:/usr/bin' }

use strict;
use warnings;
use Getopt::Long;
use Config::Tiny;
use File::Find;
use Text::Glob qw/glob_to_regex/;

my $Usage   = "$0 [--config FILE] [--dont] [--verbose] [--no-mime] FILES-AND-DIRS...\n";
my $Config  = "$ENV{HOME}/.subversion/config";
my $Dont    = 0;
my $Verbose = 0;
my $NoMIME  = 0;
GetOptions('config=s' => \$Config,
	   'dont'     => \$Dont,
	   'verbose+' => \$Verbose,
	   'no-mime'   => \$NoMIME,
) or die $Usage;

my $tiny = Config::Tiny->read($Config)
    or die "Error while reading \"$Config\": ", Config::Tiny::errstr(), "\n";

exists $tiny->{'auto-props'}
    or die "Config file $Config doesn't have an auto-props section\n";

my @autoprops;

while (my ($glob, $props) = each %{$tiny->{'auto-props'}}) {
    my $regex = glob_to_regex($glob);
    my %props;
    foreach (split /\s*;\s*/, $props) {
	my ($prop, $value) = split /=/;
	$props{$prop} = $value;
    }
    push @autoprops, [$regex => \%props];
}

sub fixprops {
    my ($file, $props) = @_;
    my $has_mime = 0;
    while (my ($prop, $value) = each %$props) {
	$has_mime ||= $prop eq 'svn:mime-type';
	$value = '*' unless defined $value;
	my $curval = `svn pg "$prop" "$file"`;
	return if $? != 0;
	my $NLs    = chomp $curval;
	if ($value ne $curval || $NLs == 0) {
	    my $cmd = "svn ps '$prop' '$value' '$file'";
	    if ($Verbose > 1) {
		warn "$File::Find::name: setting property '$prop' to '$value' ",
		    ($NLs == 0) ? "(it wasn't set)\n" : "(it was '$curval')\n";
	    }
	    print "$cmd\n" if $Verbose;
	    system($cmd) unless $Dont;
	}
    }
    $has_mime or warn "$File::Find::name: miss svn:mime-type in config file.\n";
}

sub wanted {
    unless (-f $_) {
	$File::Find::prune = 1
	    if -d _ && ($_ eq '.svn' || ! -d "$_/.svn");
	return;
    }

    foreach my $auto (@autoprops) {
	fixprops($_, $auto->[1]) if /$auto->[0]/;
    }
}

find(\&wanted, @ARGV);


__END__
=head1 NAME

svn-fix-props.pl - Fix Subversion properties on a working area.

=head1 SYNOPSIS

svn-fix-props.pl OPTIONS path...

=head1 DESCRIPTION

This script finds all files recursively under the specified paths,
checking if their Subversion properties are consistent with the ones
specified in the [auto-props] section in the Subversion configuration
file (usually $HOME/.subversion/config). Missing or wrong properties
are reset to their expected values.

The specified paths must be files or directories in Subversion working
areas.

After running this script you'll need to explicitly commit the changes
to the repository or to revert them.

=head1 OPTIONS

=over

=item --config FILE

Specify the configuration file from which to get the [auto-props]
specification. (By default: $HOME/.subversion/config.)

=item --dont

Just check but don't change anything.

=item --verbose

Say what's being done. Repeat it to get more verbosity.

=back

=head1 COPYRIGHT

Copyright 2008 CPqD.

=head1 AUTHOR

Gustavo Chaves <gustavo@cpqd.com.br>
