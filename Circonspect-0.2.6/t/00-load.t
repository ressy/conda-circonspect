#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Circonspect' );
}

diag( "Testing Circonspect $Circonspect::VERSION, Perl $], $^X" );
