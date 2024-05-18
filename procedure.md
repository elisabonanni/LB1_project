# Procedure to build an Hidden Markov Model for the detection of the Kunitz Domain starting from structural data
## 1. COLLECTION OF REPRESENTATIVE DATASETS
### 1.1 Structure Selection
To select the protein to train the Hidden Markov Model run an Advance Search on PDB:

**QUERY**: ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Data Collection Resolution <=3 AND Polymer Entity Sequence Length = [ 50 - 80 ] AND Polymer Entity Mutation Count = 0

**RESULT**: 131 structures

Grouped by 70% of similarity:

**RESULT**: 21 structures

Group the results in a Customize Tabular Report with:

Entry ID	Resolution (Å)	Sequence	Auth Asym ID	Entity ID Entry Id (Polymer Entity Identifiers)

Download it as a csv file (_rcsb_pdb_custom_report.csv_).

From the csv file retrieve Entry ID, Auth Asym ID(chain), Sequence for each result and redirect them in a new file (_pdb_seq.fasta_).
:
```
cat rcsb_pdb_custom_report.csv| tr -d '"'|tail -n +3| awk -F "," '{if ($1!=""){print ">"$6"_"$4"\n"$3}}' > pdb_seq.fasta
```

From the csv file retrieve the Entry ID and redirect them in a new file (_pdb.ids_):

```
cat rcsb_pdb_custom_report.csv | tr -d '"' | tail -n +3 | awk -F "," '{if ($1!=""){print $6":" $4}}' >pdb.ids
```

### 1.2 Positive Set
To select a positive set of protein annotated with Kunitz domain that will be used during the model test, perform an Advance Search on UniProtKB/Swiss-Prot.

**QUERY**: :(reviewed:true) AND (xref:pfam-PF00014)

**RESULT**: 391

Download the result and unzip it (_pdb_seq.fasta_).

```
zcal -f all_bpti.fasta
```

Create a file with all the ids from _all_bpti.fasta_:
```
grep ">" all_bpti.fasta| cut -d " " -f 1| tr -d ">"| sort >all.ids
```

The retrieved sequences comprise also the ones equal or similar that will be used to build the model (see _1.1 Structure Selection_). These must be removed with a BLAST search.

Create a blast database with the 21 sequences of the protein that will be used for the multiple structure alignment.

```
makeblastdb -in pdb_seq.fasta -dbtype prot
```

Run the blast using as query the file with all the positive sequences (_all_bpti.fasta_) and as database the file just created (_pdb_seq.fasta_):

```
blastp -query all_bpti.fasta -db pdb_seq.fasta -out all_bpti.blast -outfmt 7
```

The results with 100% of similarity are checked in the output file (_all_bpti.blast_).
```
grep -v "^#" all_bpti.blast|awk '{if($3 == 100 && $5==0 && $6==0) print $0}'|sort -nk 4|cut -f 1|sort -u|wc

# 20      20     414
```

The sequences that match at the 100% are 20, but for the model we have selected 21. It's better to be less strict with the similarity threshold to be sure to not bias the result. A threshold of 95% over a alignment long at least 50 residues it's tested:

```
 grep -v "^#" all_bpti.blast|awk '{if ($3>95 && $4>50) {print $0}}'|sort -nk 4 |cut -f 1 |sort -u |wc

# 33      33     685
```

The resulting sequences to remove are now 33, their ids can be grepped and stored in a new file (_tobe_removed.ids_):

```
grep -v "^#" all_bpti.blast|awk '{if($3 >95 && $4>50) print $0}'|cut -f 1|sort -u >tobe_removed.ids
```

Check how many ids are not in common with _all.ids_:

```
comm -23 all.ids tobe_removed.ids|wc

# 358     358    7679
```

Store them in a new file (_pos_selected.ids_):

```
comm -23 all.ids tobe_removed.ids| cut -d "|" -f 2  > pos_selected.ids
```

Retrieve the sequences of these selected ids from the fasta (_all_bpti.fasta_) file using the python script _select_fasta.py_:

```
python3 select_fasta.py all_bpti.fasta pos_selected.ids >bpti_pos_selected.fasta
```

### 1.3 NEGATIVE SET
The negative set of proteins without the Kunitz domain can be directly downloaded as the result of an Advance Search in UniProtKB (_negatives.fasta_):

**QUERY**: :(reviewed:true) NOT (xref:pfam-PF00014)

**RESULT**: 570,891

## 2. MULTIPLE STRUCTURE ALIGNMENT
A multiple structure alignment with the proteins selected _pdb.ids_ can be run using PDEeFold. The result is downloaded in _kunitz_3d.aln_ .

The sequences in this file needs to be written as a single line to make the following steps easier (_clean_kunitz_3d.aln_):
```
awk '{if(substr($0,1,1)==">") {print "\n"$1} else {printf "%s",$1}}' kunitz_3d.aln > clean_kunitz_3d.aln
```

Looking at this file we can see that in the initial position of the alignment there are many gaps:
```
head clean_kunitz_3d.aln| less

>PDB:1aap:A
--------------------vrevcseqaetgpcrAmISRWYFDVTEGKCAPFFYGGCGGNRNNFDTEEYCMAVCg---
>PDB:1bun:B
------------------rkrhpdcdkppdtkicqTvVRAFYYKPSAKRCVQFRYGGCNGNGNHFKSDHLCRCECleyr
>PDB:1dtx:A
------------------eprrklcilhrnpgrcyDkIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCig--
>PDB:1fak:I
--------------------apdfcleppydgpcrAlHLRYFYNAKAGLCQTFYYGGCLAKRNNFESAEDCMRTC----
```

These positions needs to be removed to avoid not having statistics in the HMM, since at least some residues are necessary to define the matching states.

In order to to that, a preliminary HMM can be built using hmmbuild program from HMMER package to understand in which position the alignment without gaps starts.

```
hmmbuild clean_kunitz_3d.hmm clean_kunitz_3d.aln

less clean_kunitz_3d.hmm

HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.71574  2.61998  3.05842  2.64921  2.84558  2.72052  3.80314  3.43106  2.58955  3.03621  4.06750  2.72625  3.50132  3.03093  2.84545  2.76158  2.91642  3.19592  4.56318  2.92048
          2.68658  4.42265  2.77523  2.73014  3.46394  2.40540  3.72386  3.29282  2.67753  2.69361  4.24673  2.90387  2.73727  3.18127  2.89813  2.37914  2.77487  2.98558  4.58398  3.61543
          0.83264  1.05753  1.52422  1.25841  0.33422  0.00000        *
      1   2.54233  5.20231  2.40772  2.08178  4.54239  3.47189  3.66300  4.02474  2.21299  3.12374  4.24493  2.70674  3.31980  2.32167  2.54099  2.66684  2.93462  3.59040  5.63699  4.22902     20 e - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.01681  4.48996  5.21231  0.61958  0.77255  0.52296  0.89835

```

Looking at the output file _clean_kunitz_3d_ it can be determined that the positions to consider start from the 20 residue, so the Multiple Structure Alignment is cleaned considering just the 20th position and the 59 after it (corresponding to the length of the alignment):

```
awk '{if(substr($0,1,1)==">") {print $0} else {print toupper (substr($0,20,59))}}' clean_kunitz_3d.aln > cut_kunitz_3d.aln
```

Check if the new alignment file (_cut_kunitz_3d.aln_) looks better:
```
head cut_kunitz_3d.aln|less

>PDB:1aap:A
-VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFDTEEYCMAVCG--
>PDB:1bun:B
KRHPDCDKPPDTKICQTVVRAFYYKPSAKRCVQFRYGGCNGNGNHFKSDHLCRCECLEY
>PDB:1dtx:A
PRRKLCILHRNPGRCYDKIPAFYYNQKKKQCERFDWSGCGGNSNRFKTIEECRRTCIG-
>PDB:1fak:I
-APDFCLEPPYDGPCRALHLRYFYNAKAGLCQTFYYGGCLAKRNNFESAEDCMRTC---
```

## 3. GENERATION OF THE HMM
From the cleaned and cut alignment (_cut_kunitz_3d.hmm_) can be finally build the HMM using hmmbuild from HMMER package:

```
hmmbuild cut_kunitz_3d.hmm cut_kunitz_3d.aln
```

## 4. MODEL TESTING
A cross-validation with a basic 2 folds can be used to test the performance. The overall positive (see _1.2 Positive Set_) representatives are split firstly in 2 groups (with a proportion of 7:3), then the larger group is divided again by 2. These last 2 sets are used to determine the best e-value threshold, while the other set will be used as test set for the evaluation of the performance. The same is done with the negatives.

### 4.1 Obtaining the sets
Sort randomly the positive selected ids (_pos_selected.ids_):
```
sort -R pos_selected.ids > r1_selected.ids
```

Split them in 2 with proportion 7:3 :
```
head -n 250 r1_selected.ids > pos_val.ids
tail -n 108 r1_selected.ids > pos_test.ids
```

Divide the larger group (_pos_test.ids_) in 2 halves:
```
head -n 125 pos_val.ids > pos_1.ids
tail -n 125 pos_val.ids > pos_2.ids
```

Check that the sets don't overlap:
```
comm -12 <(sort pos_val.ids) <(sort pos_test.ids) |wc
# 0 0 0

comm -12 <(sort pos_1.ids) <(sort pos_2.ids) |wc
# 0 0 0
```

Select the fasta sequences for all the sets of ids obtained and put them in new files _.fasta_ using the python script _select_fasta.py_:

```
python3 select_fasta.py bpti_pos_selected.fasta pos_1.ids 1 > pos_1.fasta
```

Check to have the correct number inside the set (125):
```
grep '>' pos_1.fasta|wc
#     125     125    1004
```

```
python3 select_fasta.py bpti_pos_selected.fasta pos_2.ids 1 > pos_2.fasta

grep '>' pos_2.fasta|wc
#     125     125    1016

python3 select_fasta.py bpti_pos_selected.fasta pos_test.ids 1 > pos_test.fasta

grep '>' pos_test.fasta|wc
#     108     108     872
```

The same kind of divisions can be done on the negatives, taking into consideration that the file is zipped (_negatives.fasta_):

```
zcat -f negatives.fasta | grep "^>"| cut -d "|" -f 2 >negatives.ids
```

Check the number of ids:
```
wc negatives.ids
#   570891  570891 4014789
```

Do the random shuffling:
```
sort -R negatives.ids > r1_negatives.ids
```

Split in 2 with proportion 7:3 :
```
head -n 400000 r1_negatives.ids > neg_val.ids
tail -n 170891 r1_negatives.ids > neg_test.ids
```

Divide the larger in 2:
```
head -n 200000 neg_val.ids > neg_1.ids
tail -n 200000 neg_val.ids > neg_2.ids
```

Check to not have overlapping:
```
comm -12 <(sort neg_1.ids) <(sort neg_2.ids) |wc
#0 0 0
comm -12 <(sort neg_val.ids) <(sort neg_test.ids) |wc
# 0 0 0
```

Retrieve the sequences from the ids using the python script _select_fasta.py_:

```
python3 select_fasta.py <(zcat -f negatives.fasta) neg_1.ids > neg_1.fasta

grep '>' neg_1.fasta|wc
# 200000  200000 1606804

python3 select_fasta.py <(zcat -f negatives.fasta) neg_2.ids > neg_2.fasta

grep '>' neg_2.fasta|wc
#200000  200000 1606480

python3 select_fasta.py <(zcat -f negatives.fasta) neg_test.ids > neg_test.fasta

grep '>' neg_test.fasta|wc
#  170891  170891 1372396
```

Perform the hmmsearch using HMMER package on the different positive and negative groups of Set1 and Set2:
```
hmmsearch -Z 1 --domZ 1 --noali --max --tblout pos_1.out cut_kunitz_3d.hmm pos_1.fasta
hmmsearch -Z 1 --domZ 1 --noali --max --tblout pos_2.out cut_kunitz_3d.hmm pos_2.fasta
hmmsearch -Z 1 --domZ 1 --noali --max --tblout neg_1.out cut_kunitz_3d.hmm neg_1.fasta
hmmsearch -Z 1 --domZ 1 --noali --max --tblout neg_2.out cut_kunitz_3d.hmm neg_2.fasta
```

Prepare the prediction files with: id, e-value and label (0 for the negatives and 1 for the positives):

```
grep -v "^#" pos_1.out | awk '{print $1"\t"$8"\t"1}' > set_1.txt
grep -v "^#" pos_2.out | awk '{print $1"\t"$8"\t"1}' > set_2.txt
```

Check:
```
wc set_1.txt
#     125     375    2097 set_1.txt

wc set_2.txt
#     125     375    2111 set_2.txt
```

For the negatives must be considered that not all of them will be present, because some will have probably a too high e-value. So, the also the missing ones must be retrieved and added to the prediction file:

Initialise the prediction file with the negatives retrieved by the hmmsearch:
```
grep -v "^#" neg_1.out | awk '{print $1"\t"$8"\t"0}' |sort -grk 2 > tmp_neg_1.txt
grep -v "^#" neg_2.out | awk '{print t$1"\t"$8"\t"0}' |sort -grk 2 > tmp_neg_2.txt
```

Check:
```
wc tmp_neg_1.txt
#   91318  273954 1182352 tmp_neg_1.txt

wc tmp_neg_2.txt
#   91136  273408 1179296 tmp_neg_2.txt
```

In this way it's confirmed that some of them are missing (since they should be 200000 for each set). Add the missing ones:
```
comm -23 <(sort neg_1.ids) <(cut -f 1 tmp_neg_1.txt | sort) | awk '{print $1"\t10\t0"}' >> tmp_neg_1.txt
comm -23 <(sort neg_2.ids) <(cut -f 1 tmp_neg_2.txt | sort) | awk '{print $1"\t10\t0"}' >> tmp_neg_2.txt
```

Check:
```
wc tmp_neg_1.txt
#  200000  600000 2489868 tmp_neg_1.txt

wc tmp_neg_2.txt
#  200000  600000 2488796 tmp_neg_2.txt
```

Put  together positive and negative of the corresponding sets:
```
cat tmp_neg_1.txt >> set_1.txt
cat tmp_neg_2.txt >> set_2.txt
```

Check:
```
wc set_1.txt
#  200125  600375 2491965 set_1.txt

#wc set_2.txt
#  200125  600375 2490907 set_2.txt
```

### 4.2 Optimize on Set 1 and Set 2
Use the python script _performance.py_ to find the optimal e-value threshold on a range between 1e-01 and 1e-15, first on Set 1 then Set 2 (_set_1.txt_ and _set_2.txt_):

```
for i in `seq 1 15`; do python3 performance.py set_1.txt 1e-$i; done > set_1.res
for i in `seq 1 15`; do python3 performance.py set_2.txt 1e-$i; done > set_2.res
```

Check the thresolds obtained from the output files (_set_1.res_ and _set_2.res_):

```
less set_1.res

TH: 0.1 Q2: 0.9585559025608994 MCC: 0.11929647376120928 F1: 0.029260299625468167 TP= 125.0 TN= 191706.0 FP= 8294.0 FN= 0.0
TH: 0.01 Q2: 0.9965221736414741 MCC: 0.389517018777324 F1: 0.2642706131078224 TP= 125.0 TN= 199304.0 FP= 696.0 FN= 0.0
TH: 0.001 Q2: 0.9996702061211743 MCC: 0.8088475092337969 F1: 0.7911392405063291 TP= 125.0 TN= 199934.0 FP= 66.0 FN= 0.0
TH: 0.0001 Q2: 0.999950031230481 MCC: 0.9622263920874493 F1: 0.9615384615384615 TP= 125.0 TN= 199990.0 FP= 10.0 FN= 0.0
TH: 1e-05 Q2: 0.9999750156152405 MCC: 0.9805684183558652 F1: 0.9803921568627451 TP= 125.0 TN= 199995.0 FP= 5.0 FN= 0.0
TH: 1e-06 Q2: 0.9999850093691443 MCC: 0.9882043571865587 F1: 0.9881422924901185 TP= 125.0 TN= 199997.0 FP= 3.0 FN= 0.0
TH: 1e-07 Q2: 0.9999950031230481 MCC: 0.9960213510492794 F1: 0.9960159362549801 TP= 125.0 TN= 199999.0 FP= 1.0 FN= 0.0
TH: 1e-08 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-09 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-10 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-11 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-12 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-13 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-14 Q2: 0.9999800124921924 MCC: 0.9838600715483844 F1: 0.983739837398374 TP= 121.0 TN= 200000.0 FP= 0.0 FN= 4.0
TH: 1e-15 Q2: 0.9999800124921924 MCC: 0.9838600715483844 F1: 0.983739837398374 TP= 121.0 TN= 200000.0 FP= 0.0 FN= 4.0'''



less set_2.res

TH: 0.1 Q2: 0.9582860712054966 MCC: 0.11889896751188557 F1: 0.02907652942544778 TP= 125.0 TN= 191652.0 FP= 8348.0 FN= 0.0
TH: 0.01 Q2: 0.9965171767645222 MCC: 0.3892790375807008 F1: 0.26399155227032733 TP= 125.0 TN= 199303.0 FP= 697.0 FN= 0.0
TH: 0.001 Q2: 0.9996851967520299 MCC: 0.8152816541121938 F1: 0.7987220447284346 TP= 125.0 TN= 199937.0 FP= 63.0 FN= 0.0
TH: 0.0001 Q2: 0.9999200499687695 MCC: 0.9415168085112553 F1: 0.9398496240601504 TP= 125.0 TN= 199984.0 FP= 16.0 FN= 0.0
TH: 1e-05 Q2: 0.9999850093691443 MCC: 0.9882043571865587 F1: 0.9881422924901185 TP= 125.0 TN= 199997.0 FP= 3.0 FN= 0.0
TH: 1e-06 Q2: 0.9999900062460962 MCC: 0.9920897771795918 F1: 0.9920634920634921 TP= 125.0 TN= 199998.0 FP= 2.0 FN= 0.0
TH: 1e-07 Q2: 0.9999900062460962 MCC: 0.991995 F1: 0.992 TP= 124.0 TN= 199999.0 FP= 1.0 FN= 1.0
TH: 1e-08 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-09 Q2: 0.9999950031230481 MCC: 0.9959894778685164 F1: 0.9959839357429718 TP= 124.0 TN= 200000.0 FP= 0.0 FN= 1.0
TH: 1e-10 Q2: 0.9999900062460962 MCC: 0.9919627816094709 F1: 0.9919354838709677 TP= 123.0 TN= 200000.0 FP= 0.0 FN= 2.0
TH: 1e-11 Q2: 0.9999850093691443 MCC: 0.9879197134482117 F1: 0.9878542510121457 TP= 122.0 TN= 200000.0 FP= 0.0 FN= 3.0
TH: 1e-12 Q2: 0.9999850093691443 MCC: 0.9879197134482117 F1: 0.9878542510121457 TP= 122.0 TN= 200000.0 FP= 0.0 FN= 3.0
TH: 1e-13 Q2: 0.9999750156152405 MCC: 0.9797836498941922 F1: 0.9795918367346939 TP= 120.0 TN= 200000.0 FP= 0.0 FN= 5.0
TH: 1e-14 Q2: 0.9999600249843847 MCC: 0.9674515809576932 F1: 0.9669421487603307 TP= 117.0 TN= 200000.0 FP= 0.0 FN= 8.0
TH: 1e-15 Q2: 0.9999600249843847 MCC: 0.9674515809576932 F1: 0.9669421487603307 TP= 117.0 TN= 200000.0 FP= 0.0 FN= 8.0'''
```

### 4.3 Performance Evaluation on Test Set
Take the E-Value of the higher MCC for both the outputs and do the mean. Use this e-value to test the performance on the Test Set.
In this case the E-Value to use is: (1e-07+1e-09)/2 = 1e-08

Do the same procedure as before to prepare the prediction files for the test set:

```
hmmsearch -Z 1 --domZ 1 --noali --max --tblout pos_test.out cut_kunitz_3d.hmm pos_test.fasta
hmmsearch -Z 1 --domZ 1 --noali --max --tblout neg_test.out cut_kunitz_3d.hmm neg_test.fasta

grep -v "^#" pos_test.out | awk '{print $1"\t"$8"\t"1}' > test.txt
grep -v "^#" neg_test.out | awk '{print $1"\t"$8"\t"0}' |sort -grk 2 > tmp_neg_test.txt
comm -23 <(sort neg_test.ids) <(cut -f 1 tmp_neg_test.txt | sort) | awk '{print $1"\t10\t0"}' >> tmp_neg_test.txt
cat tmp_neg_test.txt >> test.txt
```

Do the final evaluation of the performance using the python script _performance.py_ and using the e-value threshold found (1e-08):
```
python3 performance.py test.txt 1e-8 > test.res

TH: 1e-08 Q2: 0.9999941520125849 MCC: 0.9953566914798339 F1: 0.9953488372093023 TP= 107.0 TN= 170891.0 FP= 0.0 FN= 1.0
```

Check which ids correspond to the False Negative (FN) taking the result with the lowest e-value:

```
awk '{if ($3==1 && $2>1e-08) print $0}' test.txt|less
# Q8WPG5  1.6e-08 1
```
