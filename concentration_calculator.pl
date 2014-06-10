#!/usr/bin/perl

use Getopt::Long;

&Getopt::Long::Configure('pass_through', 'no_autoabbrev');
&Getopt::Long::GetOptions('help|h' => \$help,
	'absorbance|a=f' => \$absorbance,
	'extinction|e=f' => \$extinction,
	'dilution|d:f' => \$dilution,
	'pathLength|p:f' => \$pathLength,
);

if($help == 1){
	print "Calculate extinction coefficient from following parameters \n
	absorbance|a (required)= absorbance at 280 \n
	extinction|e (required)= extinction coefficient at 280 \n
	dilution|d (defaults to 1) = dilution factor \n
	path length|p (defaults to 1) = path length\n";
	exit(0);
}

if(! $absorbance){
	print "absorbance is required, use --absorbance or -a\n";
	exit(1);
}
if(! $extinction){
	print "extinction is required, use --extinction or -e\n";
	exit(1);
}
if(! $dilution){$dilution = 1;}
if(! $pathLength){$pathLength = 1;}

$concentration = (($absorbance*$dilution)/($pathLength*$extinction))*1000000;
print "Concentration = " . $concentration . "uM\n";
