#!/usr/bin/perl 

my $cmd = 'strace -e trace=mmap,munmap,brk ';
for my $arg (@ARGV) {
    $arg =~ s/'/'\\''/g;
    $cmd .= " '$arg'";
}
$cmd .= ' 2>&1 >/dev/null';
open( PIPE, "$cmd|" ) or die "Cannot execute command \"$cmd\"\n";


my $lastSize = 0;
my $currentSize = 0;
my $maxSize = 0;
my %maps;
my $topOfData = undef;
my ($addr, $length, $prot, $flags, $fd, $pgoffset);
my $newTop;
my $error = '';

while ( <PIPE> ) {
    if ( /^mmap\((.*)\) = (\w+)/ ) {
	@params = split( /, ?/, $1 );
	($addr, $length, $prot, $flags, $fd, $pgoffset) = @params;
	if ( $addr eq 'NULL' && $fd == -1 ) {
	    $maps{$2} = $length;
	    $currentSize += $length;
	}
    } elsif ( /^munmap\((\w+),/ ) {
	if ( defined( $maps{$1} )  ) {
	    $currentSize -= $maps{$1};
	    undef $maps{$1};
	}
    } elsif ( /^brk\((\w+)\)\s*= (\w+)/ ) {
	$newTop = hex( $2 );
	if ( hex( $1 ) == 0 or !defined( $topOfData ) ) {
	    $topOfData = $newTop;
	} else {
	    $currentSize += $newTop - $topOfData;
	    $topOfData = $newTop;
	}
    } else {
	print $_;
	#$error .= $_;
    }
    #if ( int( ( $currentSize - $lastSize ) / 1048576 ) != 0 ) {
#	printf( "%d\n", $currentSize / 1048576 );
#	$lastSize = $currentSize;
#    }
    if ($currentSize > $maxSize) {
	$maxSize = $currentSize;
    }

}
close( PIPE );

printf("Maximum memory used: %dM\n",$maxSize / 1048576);

if ( $? ) {
    print $error;
}

