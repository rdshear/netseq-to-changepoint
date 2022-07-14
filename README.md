# netseq-to-changepoint

### CAUTION: EXPERIMENTAL AND INCOMPLETE

A WDL workflow that estimates the location(s) of abrupt changes in RNAPII occupancy during transcriptional elongation.

This project is based on the computational pipeline described in [^1].


|          | Description                                              | Format   |
|----------|----------------------------------------------------------|----------|
| Assay    | NET-seq                                                  |          |
| Organism | *S. cerevisiae*                                          |          |
| Input    | Occupancy counts <br> seperate files for + and - strands | bedGraph |
|          | Genes (or other regions) to process                      | gff      |
| Output   | Location(s) of changepoints                              | gff3     |


[^1]: Shear, Robert D. 2020. [Inferring High Resolution Transcription Elongation Dynamics from Native Elongating Transcript Sequencing (NET-seq)](https://nrs.harvard.edu/URN-3:HUL.INSTREPOS:37365034). Master's thesis, Harvard Extension School. 
