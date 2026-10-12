# AAV9 epitope homology with common human viruses

Published AAV9 capsid epitopes stay distinct from the proteins of 50 common human viruses. Across 26 mapped peptides, the best ungapped match in those proteomes is 6 of 8 or 7 of 10 identical residues (70–75%), always with two or three substitutions. No epitope is identical to a viral peptide of length 7 or more, and none reaches 80% identity. A composition-matched shuffle of the capsid produces exact 6-mer and 7-mer overlaps at the same rate seen here, so those short matches are the background rate.

The related human parvoviruses are a separate case. Human bocavirus 1 VP1 is 36% identical to AAV9 VP1 over 610 aligned residues, and parvovirus B19 VP1 is 31% identical over 588 residues. That is family-level conservation of the capsid fold. At the homologous positions, the mapped AAV9 antibody loops and HLA peptides are not preserved: the best homologous epitope identity is 6/9 (67%) for `IPQYGYLTL` against B19. The peptides that define AAV9 serology remain AAV9-specific.

AAV2, included only as a positive control, matches AAV9 VP1 at 82% identity and reproduces the known cross-serotype T-cell peptides exactly, including `IPQYGYLTL` versus the AAV2 epitope `VPQYGYLTL`. Natural AAV infection, not the 50 viruses below, is the sequence basis for pre-existing anti-AAV immunity.

Sequence identity is a screen. It does not measure whether a T-cell receptor or antibody would cross-react, and it cannot see discontinuous antibody surfaces.

## Epitopes

The query is AAV9 capsid VP1 (UniProt Q6JC40, 736 amino acids), confirmed against UniProt before the search. The 26 peptides are published linear epitopes or the linear residue stretches of mapped conformational epitopes:

- Naturally presented HLA class I peptides of AAV9, including the previously validated ligands `IPQYGYLTL` (HLA-B*07:02), `LIDQYLYYL` (HLA-A*02:01), `QPAKKRLNF`, and `SQAVGRSSF` (Brito-Sierra et al., Front. Immunol. 2023, PMID 37662941).
- `FQFSYEFENV` (HLA-A*02:01) and `SVAGPSNMAVQGR` (HLA-A*68:01), the dominant responses in the 2025 immunopeptidome deconvolution (doi:10.3389/fimmu.2025.1641289).
- Antibody contact windows on the three-fold protrusions: VR-IV `TINGSGQNQQT` (residues 450–460, including the HL2370/HL2374 contact S454; doi:10.1128/JVI.01251-21) and the PAV9.1 contacts `TTVTQNNN` (491–498) and `SAQAQAQTG` (586–594; doi:10.1128/JVI.01011-18). The AAV1 loop sequences named in the escape-mutant paper (`QNQSGSAQN`, `KTKTD`, `SSTDP`) are the replacement sequences, not the AAV9 residues.
- Two overlapping mouse CD4 15-mers centered on R312 (doi:10.1038/s41467-026-69917-9). These are mouse epitopes.

Coordinates were checked against VP1. Every peptide in `data/epitopes.tsv` is an exact substring.

## The 50 viruses

Inclusion used lifetime seroprevalence where that is meaningful, and near-universal or repeated infection where antibodies wane (respiratory and enteric viruses). One reference proteome was used for each virus, usually an NCBI RefSeq genome, with every coded protein included. Segmented viruses include every segment. Dengue is represented by serotype 1. The accessions are in `data/viruses.tsv`.

| Group | Viruses |
|---|---|
| Herpesviruses | EBV, HHV-6B, HHV-7, VZV, CMV, HSV-1, HSV-2 |
| Polyomaviruses | BK, JC, Merkel cell, KI, WU |
| Parvoviruses | Parvovirus B19, human bocavirus 1 |
| Anellovirus | Torque teno virus 1 |
| Adenoviruses | Human adenovirus C, human adenovirus B1 |
| Picornaviruses | Rhinovirus A, B, and C; enterovirus A, B, and D68; parechovirus 1; hepatitis A |
| Influenza | A/H1N1pdm, A/H3N2, B |
| Paramyxoviruses and pneumoviruses | RSV A, metapneumovirus, parainfluenza 1, 2, 3, and 4 |
| Coronaviruses | 229E, OC43, NL63, HKU1, SARS-CoV-2 |
| Gastroenteritis viruses | Norovirus GI and GII, rotavirus A, astrovirus, sapovirus |
| Papillomaviruses | HPV-16, HPV-6b |
| Others | Hepatitis B, hepatitis E, dengue virus 1, molluscum contagiosum virus |

Together these proteomes contain 1,150 proteins and 570,212 amino acids. HIV-1, hepatitis C, measles, mumps, rubella, and HHV-8 were left out because their global prevalence or current incidence is lower than this set. AAV2 was scored separately as a control and is not one of the 50.

## Scoring

Thresholds were fixed before looking at the hits.

- **Exact:** identical peptide, length at least 7.
- **High:** length at least 8, at least 80% identity, at most 2 substitutions.
- **Moderate:** length at least 8, at least 70% identity, at most 3 substitutions.

Every epitope was also scored by its single best ungapped window in each proteome, whether or not it crossed a threshold. In parallel, every AAV9 6-mer and 7-mer was looked up exactly, and every AAV9 9-mer was compared at 7 or more identities. A null model shuffled AAV9 VP1 25 times with its amino-acid composition preserved (seed 20261012) and repeated the 6-mer, 7-mer, and 9-mer counts. Parvovirus capsids were aligned to AAV9 with BLOSUM62 local alignment so that epitope identity could be read at the homologous position rather than at the best chance window anywhere in the proteome.

## Best match for each epitope

The AAV2 column is the positive control, not one of the 50. Several antibody-loop peptides also differ from AAV2; those loops are serotype-specific, and the control's value is the exact recovery of the conserved T-cell peptides plus the 82% capsid alignment.

| Epitope | Best match in the 50 | Identity | Viral peptide | AAV2 |
|---|---|---|---|---|
| `DTSFGGNLGR` | HCoV-NL63 ORF1ab | 7/10 (70%) | `DLSFDKNLGR` | 10/10 |
| `TSFGGNLGR` | VZV UL24 | 6/9 (67%) | `MSFEVNLGR` | 9/9 |
| `QPAKKRLNF` | VZV large tegument protein | 6/9 (67%) | `KPAKKVLKF` | 8/9 |
| `RVITTSTRTW` | CMV UL8 | 6/10 (60%) | `IVTTTSTVTD` | 10/10 |
| `ITTSTRTW` | HPV-16 L2 | 6/8 (75%) | `ITTSTDTT` | 8/8 |
| `NNWGFRPKRLNFKLF` | EBV nuclear egress protein | 6/15 (40%) | `YHLEFVYKFLAFKLK` | 15/15 |
| `GFRPKRLNFKLFNIQ` | EBV nuclear egress protein | 7/15 (47%) | `EFVYKFLAFKLKNCN` | 15/15 |
| `TDSDYQLPYVL` | EBV UL17 | 6/11 (55%) | `TTSDGQAFLVL` | 10/11 |
| `SDYQLPYVL` | HCoV-229E ORF1ab | 6/9 (67%) | `ADYILPLVL` | 8/9 |
| `IPQYGYLTL` | B19 VP1 | 6/9 (67%) | `PPQYAYLTV` | 8/9 (`VPQYGYLTL`) |
| `SQAVGRSSF` | SARS-CoV-2 nucleocapsid | 6/9 (67%) | `SQASSRSSS` | 9/9 |
| `FQFSYEFENV` | HHV-7 large tegument protein | 6/10 (60%) | `KEFSEEFENA` | 7/10 |
| `LIDQYLYYL` | SARS-CoV-2 ORF7b | 6/9 (67%) | `LIDFYLCFL` | 9/9 |
| `TINGSGQNQQT` | VZV VP13/14 | 7/11 (64%) | `TYNGSEQNTST` | 4/11 |
| `SVAGPSNMAVQGR` | VZV DNA polymerase | 7/13 (54%) | `RRAMPSGMAGFGR` | 6/13 |
| `VAGPSNMAV` | HSV-1 ICP4 | 6/9 (67%) | `VAGGSEAAV` | 4/9 |
| `TTVTQNNN` | Adenovirus C penton | 6/8 (75%) | `TTVIQNND` | 4/8 |
| `SEFAWPGASSW` | HSV-1 UL56 | 6/11 (55%) | `SEAAQPDAGLW` | 5/11 |
| `FAWPGASSW` | EBV large tegument protein | 5/9 (56%) | `FNAAGASRW` | 4/9 |
| `RFFPLSGSLIF` | HSV-2 VP11/12 | 7/11 (64%) | `RFVALSGSLEA` | 8/11 |
| `RFFPLSGSLI` | HSV-2 VP11/12 | 7/10 (70%) | `RFVALSGSLE` | 7/10 |
| `SAQAQAQTG` | Adenovirus B 13.6K | 6/9 (67%) | `SAQAGAVQG` | 4/9 |
| `LPGMVWQDRDVY` | EBV DNA polymerase | 6/12 (50%) | `LPMPSASDRDGY` | 12/12 |
| `MVWQDRDVY` | HHV-6B | 6/9 (67%) | `MVDQDRKKY` | 9/9 |
| `TPVPADPPTAF` | HSV-1 DNA polymerase | 7/11 (64%) | `TPSPADPPGGA` | 8/11 |
| `VPADPPTAF` | EBV BHLF1 | 6/9 (67%) | `GPADPPAAA` | 6/9 |

The six moderate calls are `DTSFGGNLGR`, `ITTSTRTW`, `TTVTQNNN` (three viruses), and `RFFPLSGSLI`. None of them contains an identical run of 6 residues. `TTVTQNNN` is a short, asparagine-rich loop, and its adenovirus and HPV-16 matches (`TTVIQNND`, `TTVVQNND`, `TTVTTHNN`) share that composition. At the homologous position in B19 and bocavirus, this same loop is 1/8 and 2/8 identical.

`IPQYGYLTL` against B19 `PPQYAYLTV` is the closest biologically placed hit: it is the homologous capsid site, and it still differs at three positions.

## Whole-capsid k-mers

| Comparison | AAV9 VP1 | Shuffle mean (max) |
|---|---|---|
| Exact 7-mers shared with the 50 proteomes | 0 | 0.72 (3) |
| Exact 6-mers | 9 | 12.0 (18) |
| 9-mers with at least 7/9 identity | 55 windows, 60 pairs | 21.0 (32) windows |

Exact 6-mers and 7-mers are within the shuffle distribution. The 9-mer count is higher than every shuffle. That excess is concentrated in two parvovirus capsids: 29 of the pairs are in bocavirus VP1 and 5 are in B19 VP1. Five of the pairs reach 8/9, all of them overlapping windows in bocavirus VP1, and none of them falls inside a published epitope.

The two conserved stretches are:

```
AAV9 43-79   RGLVLPGYKYLGPGNGLDKGEPVNAADAAALEHDKAY
Bocavirus 9  RGWVLPGYRYLGPFNPLDNGEPVNNADRAAQLHDHAY
```

`VLPGYKYLG` versus `VLPGYRYLG` is 8/9. This is the VP1-unique region shared by parvoviruses. Further into the capsid core:

```
AAV9 291-299      CHFSPRDWQ
Bocavirus VP1     CHFSPQDWQ
```

B19 preserves a shorter core motif, `GYSTPW`, at AAV9 residues 276–281 (`YFGYSTPWGY` versus `IMGYSTPWRY` in the alignment). The nine exact hexapeptides elsewhere in the capsid, including `AQAQAQ` from the VR-VIII peptide `SAQAQAQTG` inside molluscum contagiosum MC006L, do not exceed the shuffle rate. `AQAQAQ` is a low-complexity glutamine-alanine run; the full 9-residue epitope does not have a moderate or high match.

## Homologous position in parvovirus capsids

BLOSUM62 local alignment, identity over aligned columns:

| Virus protein | Aligned columns | Identity | AAV9 span |
|---|---|---|---|
| Bocavirus 1 VP1 | 610 | 36.1% | 43–724 |
| Bocavirus 1 VP2 | 483 | 32.9% | 195–724 |
| B19 VP1 | 588 | 31.3% | 38–698 |
| B19 VP2 | 471 | 30.6% | 212–698 |
| AAV2 VP1 (control) | 734 | 82.0% | 1–736 |

Identity of each epitope at that homologous register is in `results/parvovirus_epitope_identity.tsv`. Antibody footprints land at 11–33% identity in B19 and bocavirus. The HLA peptides are mostly in that same range. The two partial exceptions are B19 `IPQYGYLTL`/`PPQYAYLTV` (6/9) and B19 `TDSDYQLPYVL`/`VDHEYKYPYVL` (6/11), which share the C-terminal `YPYVL`.

## Reading the result

The 50-virus panel does not contain a peptide match close enough to treat any published AAV9 epitope as shared with a common unrelated virus. Matches at 70–75% with two or three substitutions, and without a shared 6-residue core, are what the null model produces from capsid-sized proteins searched against half a million residues.

Bocavirus and B19 are homologous to AAV across the capsid, at about one residue in three. Their versions of the mapped AAV9 epitopes have diverged with the surface loops. The practical source of peptides that match AAV9 T-cell epitopes exactly is another AAV serotype: AAV2 reproduces `DTSFGGNLGR`, `RVITTSTRTW`, `LIDQYLYYL`, `SQAVGRSSF`, `LPGMVWQDRDVY`, and the mouse CD4 peptides without a substitution.

## Limits

Each virus is one reference genome, not every circulating genotype. A match would have to be created by strain variation that this reference lacks. Discontinuous antibody epitopes can resemble each other without a linear peptide match. T-cell cross-reactivity is occasionally reported for peptides that differ at one or two positions; this screen finds no AAV9 epitope in the 50 proteomes at that distance. HLA binding was not rescored, and no peptide was tested experimentally.

## Reproduce

```bash
python3 src/fetch_proteomes.py   # refreshes data/proteomes from NCBI and checks UniProt Q6JC40
python3 src/homology.py          # writes results/
```

Dependencies are in `requirements.txt` (`biopython`, `numpy`, `requests`). The proteomes already in `data/proteomes/` are the CDS translations used for the tables above.

Output:

- `results/epitope_best_hits.tsv` — best viral window for each epitope
- `results/epitope_by_virus.tsv` — best window in every virus
- `results/notable_epitope_hits.tsv` — windows that cross the exact, high, or moderate bar
- `results/exact_kmers.tsv` — shared 6-mers
- `results/nine_mer_hits.tsv` — 9-mers at 7/9 or better
- `results/domain_alignments.tsv` — BLOSUM62 capsid alignments
- `results/parvovirus_epitope_identity.tsv` — epitope identity at the homologous register
- `results/summary.json` — counts and the null model
