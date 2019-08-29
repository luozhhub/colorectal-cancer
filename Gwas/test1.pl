use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  #-port => 3337
);

my $chr = 6;  #defining the region in chromosome 6
my $start = 25_834_000;
my $end = 25_854_000;


sub getLD{
my ($chr, $start, $end, $snpID) = @_;
my $population_name = '1000GENOMES:phase_3:YRI'; #we only want LD in this population

my $slice_adaptor = $registry->get_adaptor('human', 'core', 'slice'); #get adaptor for Slice object
my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$start,$end); #get slice of the region

my $population_adaptor = $registry->get_adaptor('human', 'variation', 'population'); #get adaptor for Population object
my $population = $population_adaptor->fetch_by_name($population_name); #get population object from database

my $ldFeatureContainerAdaptor = $registry->get_adaptor('human', 'variation', 'ldfeaturecontainer'); #get adaptor for LDFeatureContainer object
# Retrieve 1000 genomes phase 3 data from VCF file
$ldFeatureContainerAdaptor->db->use_vcf(1);
my $ldFeatureContainer = $ldFeatureContainerAdaptor->fetch_by_Slice($slice,$population); #retrieve all LD values in the region

open OP, "> ./YRI/$snpID" or die "failed, $1";
foreach my $r_square (@{$ldFeatureContainer->get_all_r_square_values}){
  if ($r_square->{r2} > 0.8){ #only print high LD, where high is defined as r2 > 0.8
       print OP "High LD between variations ", $snpID, "\t", $r_square->{variation1}->variation_name, "-", $r_square->{variation2}->variation_name, "\n";
         }
         }
close(OP);
}

open(IP, "< ./new.list.GWAS.site.txt") or die "failed, $!";
my @lines = <IP>;
for (my $i=0; $i < scalar(@lines); $i++){
    my $line = $lines[$i];
    chomp($line);
    my @array = split /\t/, $line;
    my $chr = $array[0];
    my $start = int($array[1]) - 250000;
    my $end = int($array[1]) + 250000;
    my $snpID = $array[2];
    getLD($chr, $start, $end, $snpID);
}
close(IP);

