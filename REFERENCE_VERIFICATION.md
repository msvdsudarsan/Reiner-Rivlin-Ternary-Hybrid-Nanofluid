# Reference verification log

The manuscript cites 43 works, and `References.bib` contains exactly those 43
entries; nothing uncited is carried in the bibliography file. This log records
how each entry was checked, so that a reader can see which citations rest on a
direct check against the publisher's record and which do not. It is kept as a
transparency record, not because the journal requires it.

## Checked in the 3 October 2026 pass

| Entry | What was checked and changed |
|---|---|
| Sadiya & Sucharitha (2025), Discover Nano 20:168 | Publisher page (Springer) and PMC record. **Corrected:** the first author is A. B. Sadiya, as the publisher's own citation line and the PMC record give. An earlier version of the bibliography printed "Sadiya" alone, and an earlier version of this log wrongly recorded that as the correct form. Volume 20, article 168 and the DOI were confirmed. |
| Mitschka & Ulbrecht (1965), Collect. Czech. Chem. Commun. 30:2511-2526 | Journal digital archive (CCCC) gives the full title, including "Ostwald-de-Waeleschen Typs in der Umgebung rotierender Drehkegel und Scheiben", issue 8, and DOI 10.1135/cccc19652511. The title, issue and DOI are now complete, and the BibTeX key reads Mitschka1965. |
| Rogers & Lance (1960), J. Fluid Mech. 7:617-631 | Year was already printed as 1960; the BibTeX key was still Rogers1961 and is now Rogers1960. DOI added. |
| Griffiths, Stephen, Bassom & Garrett (2014), JNNFM 207:1-6 | DOI 10.1016/j.jnnfm.2014.02.004 confirmed against two institutional repository records. |
| Housiadas (2026), JNNFM 349:105609 | Was cited and present in `References.bib` but missing from earlier versions of this log, which is why they said 42. |
| Zenodo archive (this work) | Version label removed from the citation so that it cannot go stale; DOI unchanged. |

DOIs added in this pass to entries that previously had none: von Karman (1921),
Cochran (1934), Turkyilmazoglu (2012, 2014), Rogers & Lance (1960), Mitschka &
Ulbrecht (1965), Weidman et al. (2006), Harris et al. (2009), Reiner (1945),
Rivlin (1948), Shevchuk (2009, book DOI), Griffiths et al. (2014, two entries),
Abdulameer et al. (2016), More et al. (2026) and Gregory, Stuart & Walker
(1955). Apart from the entries named in the table above, these DOIs were taken
from an independent audit of the bibliography dated 3 October 2026 and were not
each resolved again by the authors in this pass; every DOI in the file should be
opened once before the final upload. The informal note attached to the
Turkyilmazoglu (2012) entry was removed from the bibliography.

## Checked against the publisher record in the earlier October 2026 pass

Karami & Stephanou (2026); Salas-Barzola et al. (2026); Evans, Palhares &
Afonso (2026); Jiang, Papageorgiou & Ding (2026); Winters et al. (2026);
Rinkens et al. (2026); More, Pashkovski, Patterson & McKinley (2026, volume not
stated by the publisher, so not given); Maklad & Poole (2021); Griffiths,
Garrett & Stephen (2014); Abdulameer, Griffiths, Alveroglu & Garrett (2016);
Gregory, Stuart & Walker (1955). These were checked against the publisher
record or its BibTeX export.

## Checked in the 18 August 2026 pass

Naganthran et al. (2020); Mandal & Pal (2025); Jamrus et al. (2025); Paul,
Patgiri & Sarma (2025); Algehyne et al. (2026, Discover Nano); Gizewu & Ibrahim
(2025); Asma, Othman & Muhammad (2019, *Mathematics*); Rauf et al. (2025); Bartwal et al.
(2025). Each was located through the publisher page or DOI resolver and its
authors, journal, volume, pages and DOI compared with the citation.

## Checked in an earlier review round

Puspanathan et al. (2024); Mandal & Pal (2024); Turkyilmazoglu & Senel (2013).

## Not individually re-checked field by field

- Foundational works: von Karman (1921), Cochran (1934), Reiner (1945), Rivlin
  (1948).
- Standard method and background references: Merkin (1986), Weidman et al.
  (2006), Harris et al. (2009), Turkyilmazoglu (2012, 2014), Tabassum &
  Mustafa (2018), Lv et al. (2021), Sahoo & Kumar (2020), Shevchuk (2009), and
  Asma, Othman, Muhammad et al. (2019, *Symmetry*).

These are established, widely cited sources, so the risk of an error is low.

## The authors' own dataset

The Zenodo entry cites the archive without a version label. After a new release
is created, the DOI should be confirmed to resolve to the intended record.

## Summary

| Status | Entries |
|---|---|
| Checked in the 3 October 2026 pass, in the table above | 5 |
| Checked in the earlier October 2026 pass | 11 |
| Checked in the 18 August 2026 pass | 9 |
| Checked in an earlier review round | 3 |
| Not individually re-checked field by field | 14 |
| Authors' own dataset | 1 |
| **Total cited, and total in `References.bib`** | **43** |

The full text of the 18 August 2026 log is kept in `DEVELOPMENT_HISTORY.md`.
