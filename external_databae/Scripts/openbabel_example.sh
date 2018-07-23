# This example go through all lowest conformers and
# transform fhiaims files to PDB format using openbabel
#
#

for aa in ala arg argH asn asp aspH cys gln glu gluH gly hisD hisE hisH ile leu lys lysH met phe pro ser thr trp tyr val
do
    for cap in dipeptides uncapped
    do
	for dir in `find $aa/$cap/. -maxdepth 1 -type d`
	do
	    obabel -i fhiaims $dir/conformer.0001.fhiaims -o pdb -O $dir/conformer.0001.pdb 
	done
    done
done
