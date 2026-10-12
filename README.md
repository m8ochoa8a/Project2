# AAV9 capsid homology with common human viruses

The full AAV9 capsid does not match the surface proteins of 50 common human viruses, except for the two human parvoviruses in the set. Bocavirus 1 VP1 is 36% identical to AAV9 VP1 over 610 residues, and parvovirus B19 VP1 is 31% identical over 588. Outside those capsids, the best 15-residue window on any virion surface protein is 9 of 15 identical residues, which is the best a composition-matched shuffle of AAV9 VP1 also reaches. No surface protein contains an exact AAV9 peptide of length 7 or more, and none of the close windows falls inside a published AAV9 epitope.

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

## Full capsid against virion surface proteins

The epitope screen asked whether 26 published peptides occur anywhere in a viral proteome. This screen asks a different question: whether any part of the AAV9 capsid resembles a protein on the outside of the virus particle. The query is the whole VP1 sequence, all 736 residues, including the VP2 and VP3 regions. Rep proteins are not on the particle and were not scanned.

A surface protein here means a protein of the infectious particle that faces the outside. For a non-enveloped virus that is the outer capsid. For an enveloped virus it is the envelope glycoproteins, including the small envelope proteins. Internal capsids of enveloped viruses, tegument, polymerases, and membrane proteins that act in the infected cell rather than on the particle were left out. The catalog is `data/surface_catalog.tsv` (179 proteins, 74,856 residues).

Choices that change what gets scanned:

- Polyomavirus surface is major capsid VP1 only. JC VP1 is NP_043511.1; the FASTA product name on that record is "hypothetical protein."
- Adenovirus surface is hexon, penton base, fiber, and protein IX.
- Picornavirus outer-capsid chains were cut from the polyprotein at the NCBI mature-peptide bounds. VP4 was not taken. Parechovirus has no VP4 cleavage, so VP0 was taken. Hepatitis A chains are the UniProt P08617 bounds on NP_041007.1 (VP2 24–245, VP3 246–491, VP1 492–764); the RefSeq feature table does not list them.
- Sapovirus VP1 is ORF1 residues 1738–2280, from the Calici_coat domain through the end of the polyprotein. The minor capsid protein from ORF2 was left out.
- Rotavirus surface is outer-capsid VP4 and VP7. Norovirus, astrovirus, and hepatitis E contribute their capsid proteins.
- Hepatitis B is the large envelope protein only. It already contains preS1, preS2, and S, so the middle and small envelope proteins would have scanned the same S sequence again.
- Papillomavirus L1 is the major capsid. L2 is included and labeled as the minor capsid protein because part of it is exposed and part is internal.
- Herpesvirus surface proteins are the virion envelope glycoproteins (gB, gC, gD, gE, gG, gH, gI, gK, gL, gM, gN, and the extras each virus actually carries, such as EBV gp350/gp42/gp150 and the CMV pentamer subunits UL130, UL131A, and UL132). HSV-1 and HSV-2 each contribute 12, VZV contributes 9. UL20, UL43, and the CMV US2–US11 immune-evasion proteins are not in the set.
- Influenza surface is HA, NA, and M2. Influenza B is HA, NA, and NB; BM2 is not annotated on B/Lee/1940. Coronavirus surface is spike, envelope, membrane, and hemagglutinin-esterase where that protein exists.
- Molluscum contagiosum proteins are named only as MC###. Fifteen ORFs match a vaccinia Copenhagen virion-membrane protein at 30% or more identity across at least half of that protein. A sixteenth, MC084L, matches vaccinia H3L over 256 of 324 residues at 29.7% (score 357) and was kept because the alignment runs the length of both proteins. Two partial overlaps that cleared 30% on a gene that already had a full-length molluscum homolog (MC042L against L5R, MC046L against G3L) were dropped. No ORF reached the cutoff for vaccinia D8L, B5R, A33R, A34R, or A56R.

Every one of the 50 viruses has at least one surface sequence. AAV2 VP1 was aligned as a positive control and was not part of the null.

The same peptide rules were used. Exact means an identical peptide of length at least 7. High means length at least 8, at least 80% identity, and at most 2 substitutions. Moderate means length at least 8, at least 70% identity, and at most 3 substitutions. Every AAV9 9-mer, 11-mer, and 15-mer was also given its single best window on each virus's surface proteins. VP1 was aligned to each surface protein with BLOSUM62 (gap open −11, extension −1). The null shuffled VP1 25 times with the same composition and seed 20261012, and scored those shuffles against this surface set only.

| Comparison | AAV9 VP1 | Shuffle mean (max) |
|---|---|---|
| Exact 7-mers on surface proteins | 0 | 0.04 (1) |
| 9-mers with at least 7/9 identity | 34 windows, 36 pairs | 2.68 (7) windows |
| Best 15-mer identity anywhere on the surface set | 12/15 | 8.48 (9) |

Thirty-two of the 34 nine-residue windows are in bocavirus or B19 capsid protein. The other two are single 7/9 matches, HSV-2 glycoprotein N `GEPVNAADA` versus `GEPPNAAGA` and HHV-6B gB `YKSNNVEFA` versus `YKSVNVEEA`. Matches of that grade are what the shuffle produces. None of the 34 windows overlaps a published epitope. The 8/9 windows are the same bocavirus stretches found in the whole-proteome screen: `VLPGYKYLG` versus `VLPGYRYLG` in the VP1-unique region, and `CHFSPRDWQ` versus `CHFSPQDWQ` in the capsid core.

The best 15-residue window in each virus:

| Identity | Viruses |
|---|---|
| 12/15 | Bocavirus 1, `RGLVLPGYKYLGPGN` versus VP1 `RGWVLPGYRYLGPFN` (residues 43–57). This is the only moderate 15-mer. It has three substitutions, so it does not meet the high bar, and it is not inside a mapped epitope. |
| 10/15 | B19, `AYFGYSTPWGYFDFN` versus VP1 `PIMGYSTPWRYLDFN` (residues 273–287). This is above every shuffle but below the moderate bar of 12/15. The shared core is `GYSTPW`. |
| 9/15 | HSV-2 only. Glycoprotein E residues `SVPPPTPPRRPPVAP` match `SVPDPQPIGEPPAAP`. The shuffle maximum is also 9/15, and the stretch is proline-rich. |
| 8/15 | EBV, HHV-6B, HHV-7, HSV-1, adenovirus C, rhinovirus B, parainfluenza 3, OC43, HPV-16, hepatitis A, molluscum contagiosum |
| 7/15 | VZV, CMV, KI, WU, torque teno virus, adenovirus B, influenza A H1N1, influenza A H3N2, influenza B, metapneumovirus, parainfluenza 4, 229E, NL63, HKU1, SARS-CoV-2, norovirus GI, norovirus GII, rotavirus, astrovirus, sapovirus, enterovirus B, enterovirus D68, parechovirus, HPV-6, hepatitis E, dengue 1 |
| 6/15 | BK, JC, Merkel cell, rhinovirus A, rhinovirus C, RSV A, parainfluenza 1, parainfluenza 2, enterovirus A, hepatitis B |

BLOSUM62 local alignment separates the same two proteins from everything else. Bocavirus VP1 scores 742 (36.1% over 610 columns, AAV9 residues 43–724). B19 VP1 scores 567 (31.3% over 588 columns, residues 38–698). The highest score on any other surface protein is 61, for VZV glycoprotein I over 41 columns. AAV2 VP1, the control, scores 3294 (82.0% over 734 columns) and contains exact 15-mers of AAV9.

## Reading the result

The 50-virus panel does not contain a peptide match close enough to treat any published AAV9 epitope as shared with a common unrelated virus. Matches at 70–75% with two or three substitutions, and without a shared 6-residue core, are what the null model produces from capsid-sized proteins searched against half a million residues.

Restricting the search to virion surface proteins, and extending the query from the 26 epitopes to all 736 capsid residues, does not add a second hit. The only surface proteins with a real alignment to AAV9 VP1 are bocavirus and B19. Their versions of the mapped epitopes have diverged with the surface loops. The practical source of peptides that match AAV9 T-cell epitopes exactly is another AAV serotype: AAV2 reproduces `DTSFGGNLGR`, `RVITTSTRTW`, `LIDQYLYYL`, `SQAVGRSSF`, `LPGMVWQDRDVY`, and the mouse CD4 peptides without a substitution.

## Limits

Each virus is one reference genome, not every circulating genotype. A match would have to be created by strain variation that this reference lacks. Discontinuous antibody epitopes can resemble each other without a linear peptide match. T-cell cross-reactivity is occasionally reported for peptides that differ at one or two positions; this screen finds no AAV9 epitope in the 50 proteomes at that distance. HLA binding was not rescored, and no peptide was tested experimentally.

The surface set is a definition, not a measurement of what is exposed on a particle. L2 and a few herpes glycoproteins are only partly exposed. Molluscum proteins were assigned by homology to vaccinia membrane proteins, so a genuine molluscum surface protein with no remaining vaccinia similarity would have been missed, and a homologous protein that is not on the molluscum particle would have been included. Influenza B BM2 is absent from the RefSeq annotation used here.

## Reproduce

```bash
python3 src/fetch_proteomes.py   # refreshes data/proteomes from NCBI and checks UniProt Q6JC40
python3 src/homology.py          # epitope and whole-proteome screen; writes results/
python3 src/surface_scan.py      # full VP1 against virion surface proteins; writes results/surface_*
```

Dependencies are in `requirements.txt` (`biopython`, `numpy`, `requests`). The proteomes already in `data/proteomes/` are the CDS translations used for the tables above. `data/features/polyproteins.ft` caches the mature-peptide coordinates, and `data/reference/vaccinia_copenhagen.faa` is the vaccinia membrane-protein reference used for molluscum.

Output:

- `results/epitope_best_hits.tsv` — best viral window for each epitope
- `results/epitope_by_virus.tsv` — best window in every virus
- `results/notable_epitope_hits.tsv` — windows that cross the exact, high, or moderate bar
- `results/exact_kmers.tsv` — shared 6-mers
- `results/nine_mer_hits.tsv` — 9-mers at 7/9 or better
- `results/domain_alignments.tsv` — BLOSUM62 capsid alignments
- `results/parvovirus_epitope_identity.tsv` — epitope identity at the homologous register
- `results/summary.json` — counts and the null model
- `data/surface_catalog.tsv` — the 179 surface proteins and why each was included
- `data/surface_proteins.faa` — their sequences, with polyprotein chains already sliced
- `results/surface_best_windows.tsv` — best AAV9 9-mer, 11-mer, and 15-mer on each virus
- `results/surface_peptide_hits.tsv` — exact 7-mers and 9-mers at 7/9 or better
- `results/surface_alignments.tsv` — BLOSUM62 local alignment of VP1 to each surface protein
- `results/surface_summary.json` — surface-set counts and the surface-only null
- `results/mcv_surface_homologs.tsv` — every molluscum ORF against its best vaccinia membrane protein
