These scripts will generate a list of core genes pulled from the Sybase prokyaryotic annotation database.

You can find the genes on a web page at:

	http://www.jcvi.org/cgi-bin/genome-properties/GenomePropDefinition.cgi?prop_acc=GenProp0799

They are associated with the gene property: 

	GenProp0799

1.) run extract.csh, to pull all the core genes.
2.) run get_nonProtein_core.csh, to subset a list of all the genes not
associated with proteins

Files will be created which contain gene ids, hmms, GO and EC accessions.


###########################################################################################################

EC's with - were augmented:

Exinuclease- looks like it.s 3.1.25.1
DNA primase- we should be ok with  EC 2.7.7.6  (Catalyzes
DNA-template-directed extension of the 3'-end of an RNA strand by one
nucleotide at a time.)
tRNA ligases- I think they are all 6.1.1.X- I don.t think when HMP was looking
for the core set, all tRNA ligases ended up as .core.. For the most part, it
would be safe to consider them core as the vast majority of organisms are
going to have them (there are going to be some exceptions but that can.t be
helped)
tRNA(Ile)-lysidine synthetase- EC 6.3.4.19 is specific for lysidine
synthetase. I think this should work for now. We could choose others if that
works better for our purposes. It will be a similar situation to the tRNA
ligases in that most will be in most organisms with some exceptions that can.t
be helped
Elongation factor 4- EC 3.6.5.n1 refers to the correct enzyme so this should
work


