#!/usr/bin/perl
# Bogotá, June 8th 2015
# Bioinformatics and Systems Biology Group
# Institute for Genetics - National University of Colombia
# 
# 
# Minimal Validation for Stoichiometry Reactions 
# FUNCTIONALITY:
# For a given set of reactions, this script identifies:
# 1) All substrates, and their rnxs, that are not present as products.
# 2) All products, and their rxns, that are never used as substrates (dead ends).
#
#
# INPUT
# A two columns text file (tab separated) with the following mandatory format:
# Column 1: Reaction ID.
# Column 2: Reaction formula.
#
# CONSIDERATIONS ON INPUT FORMAT:
# This script assumes that all your reactions have the following format:
#
#            L-Glutamate[c] <=> CO2[c] + 4-Aminobutanoate[c]
#
# Where arrows and plus signs are surrounded by a "space character".
# It is also expected that stoichiometry coefficients are surrounded by spaces,
# (nothe the "2" before the CO2[c] or the NH3[c]): 
# 
#          H2O[c] + Urea-1-carboxylate[c] <=> 2 CO2[c] + 2 NH3[c]
#
# It also expects arrows to be in the form "<=", "=>" or "<=>". Meaning that arrows
# like "==>", "<==>" or "<==" will not be parsed and will lead to errors.
#
#
# OUTPUT
# This script outputs 5 files into the current folder: 
#
# INPUT.orphanProducts
#       lists all products never used as substrates. Possible dead ends.

# INPUT.orphanProducts.rxns
# 		lists all the rxns holding orphan products. 

# INPUT.orphanSubstrates
#		lists all substrates never produced in a reaction. Possible candidates to be
#		introduced into the system by exchange rxns or by adding more internal rxns.

# INPUT.orphanSubstrates.rxns
#		lists all rxns holding orphan substrates.

# INPUT.mets
#	Non redundant list of metabolites present in all rxns in input file.
#
# Maintained by: Andrés Pinzón, ampinzonv@unal.edu.co
#
#

use strict;
use warnings;
use Data::Dumper;
use File::Temp qw/ tempfile tempdir /;
use List::MoreUtils qw(uniq);


#open file
my $filename = $ARGV[0];
open( FILE, "<", $filename ) or die "Cannot open file $!";
my @data = <FILE>;

#define some vars
my @allLefts;
my @allRights;
my $fs;
my $fp;
my $dir;
my $substrate;
my $product;

#put rxnsID and rxnsFormula into a hash
my %hash = ();
foreach my $line (@data) {
    my @column = split( /\t/, $line );
    $hash{ $column[0] } = $column[1];
} 
#print Dumper(\%hash);

#and now into their own array: http://www.tutorialspoint.com/perl/perl_hashes.htm
my @rxnsForm = values %hash;
my @rxnsID   = keys %hash;



foreach my $row (@rxnsForm){
	chomp($row);
	#Get left and right sections of equation

		if($row =~ " <=> ") {
			$fs = " <=> ";
			$dir = 3;
		}elsif($row =~ " => ") {
			$fs = " => ";
			$dir = 1;
		}else {
			$fs = " <= ";
			$dir = 2;
		}


	my @sides = split($fs,$row);
	# Take care of the following cases:
	# substrate => product
	# product <= substrate
	# ? <=> ?
	# Right or left depend on the rxn direction.
	if ($dir == 1){
		$substrate = $sides[0];
		$product = $sides[1];
	}elsif($dir == 2){
		$substrate = $sides[1];
		$product = $sides[0];
	# In this case of bidirectionality of the rxns, we will assume that both sides of the
	# reaction can be both, substrate and product. For now can't imagine a better
	# solution.
	# My guess here is that at the end the mets involved in this rxns will be cancelled
	# because they are going to be at both sides of the equation?!
	#
	}else{
		$substrate = $sides[0] ." + ". $sides[1];
		$product = $sides[1] ." + ". $sides[0];
	}	

	
	
	#Get substrates
	my @left = split(/ \+ /,$substrate);
	#Clean it up
	foreach my $met (@left){
		push @allLefts, cleanup($met)."\n";
	}	 

	#Get products
	my @right = split(/ \+ /,$product);
	foreach my $met (@right){
		 push @allRights, cleanup($met)."\n";
	} 

}

# ===========================================================================
# PART 2
# To this point we have two arrays, @allLefts that holds all the substrates and
# @allRights that holds all the products. So we can start performing the analysis
# over them.

@allLefts  = sort(@allLefts);
@allRights = sort(@allRights);


#get rid of redundacy
@allLefts  = uniq(@allLefts);
@allRights = uniq(@allRights);


#get the number of substrates and products
my $numSubstrates = @allLefts;
my $numProducts   = @allRights;

#write substrates to a tmp file
($fs,my $substratesFile) = tempfile();
open (FILE, ">> $substratesFile") || die "Can not open substrates file for writing\n";
print FILE @allLefts;

#write products to a tmp file
($fp,my $productsFile) = tempfile();
open (FILE, ">> $productsFile") || die "Can not open products file for writing\n";
print FILE @allRights;


# === IDENTIFY ORPHAN SUBSTRATES  WRITE TO FILE ===
# thanks to: http://perlmeme.org/faqs/system/system.html
# returns three columns: 1)unique in substrates 2)unique in products 3)common to both.
# Here we supress columns 2 and three, so only get unique substrates.
# Note that we could also take those that are only present as products (col 2), and
# could be candidates to be dead-ends, since they are produced but never consumed.
my $uniqSubstrates = qx(comm -23 $substratesFile $productsFile);
my $uniqSubstratesFile = writeFile($filename, $uniqSubstrates );
qx(cp $uniqSubstratesFile ./$filename."orphanSubstrates");


# ====  IDENTIFY ORPHAN PRODUCTS WRITE TO FILE ====
my $uniqProducts = qx(comm -13 $substratesFile $productsFile);
my $uniqProductsFile = writeFile($filename, $uniqProducts );
qx(cp $uniqProductsFile ./$filename."orphanProducts");


# ====  IDENTIFY RXNS CONTAINING ORPHAN SUBSTRATES WRITE TO FILE ====
my @orphanSubstrate = split(/\n/,$uniqSubstrates);
my $k;
my $v;
my $orphanSubstratesRxns;

foreach my $i (@orphanSubstrate){
	chomp($i);
	while (($k, $v) = each %hash){
		chomp($v);
		
		# Quite confusing scaping the whole pattern in "=~" : 
		# http://stackoverflow.com/questions/2156731/how-do-i-escape-special-chars-in-a-string-i-interpolate-into-a-perl-regex
		# print "checking if ".$i." is contained in: ".$v."\n";
		if ($v =~ m/\Q$i\E/){
			$orphanSubstratesRxns .= $k."\t".$v ."\t".$i."\n";
		}
	}
			
}

my $orphanSubstratesRxnsFile = writeFile($filename, $orphanSubstratesRxns);
qx(cp $orphanSubstratesRxnsFile ./$filename."orphanSubstrates.rxns");


# ====  GET ALL RXNS CONTAINING ORPHAN PRODUCTS ====
my @orphanProduct = split(/\n/,$uniqProducts);
my $l;
my $w;
my $orphanProductRxns;

foreach my $j (@orphanProduct){
	chomp($j);
	while (($l, $w) = each %hash){
		chomp($w);
		
		#Quite confusing scaping the whole pattern: 
		# http://stackoverflow.com/questions/2156731/how-do-i-escape-special-chars-in-a-string-i-interpolate-into-a-perl-regex
		#print "checking if ".$i." is contained in: ".$v."\n";
		if ($w =~ m/\Q$j\E/){
			$orphanProductRxns .= $l."\t".$w ."\t".$j."\n";
		}
	}
			
}

my $orphanProductsRxnsFile = writeFile($filename, $orphanProductRxns);
qx(cp $orphanProductsRxnsFile ./$filename."orphanProducts.rxns");


# === GET A NONE REDUNDANT LIST OF METABOLITES WRITE TO FILE ===
my $metsu = qx(cat $substratesFile $productsFile| sort | uniq -u);
my $metsd = qx(cat $substratesFile $productsFile| sort | uniq -d);
my $finalMets = $metsu.$metsd;

my $finalMetsFile = writeFile($filename, $finalMets );
qx(cp $finalMetsFile ./$filename."mets");

#USER OUTPUT
my $nums = qx(wc -l $filename*);
print "Summary of analysis.\nNumber of items found or analyzed / file name\n";
print "$nums";



# === ROUTINES  ====
# This routine removes some characters from each metabolite, such
# as the compartiment strings: [c],[b], etc as well as any stoichiometric
# character.

sub cleanup {
	my ($met) =@_;
	chomp($met);
	#Remove stocihiometry coefficient
	#NOTE: this pattern is very fragile to strings that
	#do not follow strictly the rule of number+space
	$met =~s/^\d //;
	return $met;	
	#print $met."\n";
}


sub writeFile{
	my($filename, $content) = @_;
	(my $fhandle, $filename) = tempfile();
	open (FILE, ">> $filename") || die "Can not open file: $filename for writing\n";
	print FILE $content;
	return $filename;
}


exit 0;
