MUMdex python module

TODO:
1. Test all python scripts.
2. Create small size project that can be used for regression testing.
3. Write documentation.
4. Write conda recipes.
5. Modify setup, so that it will automatically recompile all c++ executables.
6. Put the repository to github.
7. Put so_tumor_mumdex pipeline to github.
8. Write new python scripts for various analysis of mums and bridges data.
9. To do 1 we need partially recreate Peter's mums repository:
   /mnt/wigclust20/data/safe/paa/analysis/mums
   the is about 0.6T.

Some missing headers:
pop2txt:
chr1 pos1 high1 chr2 pos2 high2 orientation_char(=|i|o) invariant offset n_people n_bridges median_bridges max_bridges

bridges2txt:
chr pos high chr2 pos2 high2 it inv ioff al bl aml bml bc amc bmc

population_bridges:

family parent kids nFam nSam nInFam chr chrA posA highA chrB posB highB type ori invariant offset bridge supA supB mCA mSA mCB mSB mapA mapB minParCov bridges anchorsA anchorsB coveragesA coveragesB pNP pNB pMed pMax tBC nBC nS