# Hidden species diversity within the Mediterranean cone snail, *Lautoconus ventricosus*
## What was this study about?
I am not sure if this is entirely true, but I have heard several times the Mediterranean cone snail, *L. ventricosus*, holds the records as the mollusk species with more synonym names. Regardless if that little fact is true or not, a quick look at the [WoRMS database](https://www.marinespecies.org/aphia.php?p=taxdetails&id=428401) is enough to find out that almost 100 names have been proposed for this species.

In this study, we have used *cox1* barcodes, complete mitochondrial genomes and low-coverage nuclear genomes to delve into the diversity of cone snails in the Mediterranean Sea. We propose that *L. ventricosus* is actually a complex of cryptic species, with at least three different species hiding within the same name. Some support for a fourth species is found in the data, but the signal is not strong enough and more specimens should be included before making a decision.

## Repository description
This repository was originally created to share the R scripts used to draw the main figures of this manuscript. Its intended use, however, has been expanded to provide further details to replicate the main analyses. Due to the size of the files only a handful of the input data will be provided, while all other will only be described.

### Molecular data can be downloaded from the NCBI:
<ul>
    <li><strong>Barcodes:</strong> <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=ON951339:ON951583[accn]">new sequences</a>, plus <a href="https://www.ncbi.nlm.nih.gov/nuccore/AY588229">AY588229</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KX263251">KX263251</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491607">MF491607</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KJ550006.1">KJ550006</a>; and outgroups: <a href="https://www.ncbi.nlm.nih.gov/nuccore/NC_013243">NC_013243</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801847">KY801847</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491520">MF491520</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491522">MF491522</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491523">MF491523</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491534">MF491534</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491540">MF491540</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491549">MF491549</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491565">MF491565</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY864972">KY864972</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801863">KY801863</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801862">KY801862</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801859">KY801859</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801856">KY801856</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801849">KY801849</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY864973">KY864973</a>.

  </li>
    <li><strong>Mitochondrial genomes:</strong> <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=ON968966%3AON968984%5Baccn%5D">these</a> and <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=ON975007%3AON975008%5Baccn%5D">these</a> new mitogenomes, plus <a href="https://www.ncbi.nlm.nih.gov/nuccore/KX263251">KX263251</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491607">MF491607</a>; and outgroups: <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801847">KY801847</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY801863">KY801863</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY864972">KY864972</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/KY864973">KY864973</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491520">MF491520</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491522">MF491522</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491523">MF491523</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491534">MF491534</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491540">MF491540</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/MF491549">MF491549</a>.</li>
    <li><strong>SRA raw reads:</strong> <a href="https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA856832">new low coverage genomes</a>. Please, note that because of an error the R2 reads of the sample SZN064 were never sent to us. Hence, this sample was uploaded to GenBank as single-end sequencing. In practice, if you really want to replicate our results you can duplicate this file with a new name and just use it as if it was the reverse reads. I apologize for this inconvenience.</li>
</ul>

## Citation
<ul>
  <li>Abalde, S., Crocetta, F., Tenorio, M.J., D'Aniello, S., Fassio, G., Rodriguez-Flores, P.C., Uribe, J.E., Afonso, C.M.L., Oliverio, M. and Zardoya, R. (early release). Hidden species diversity within the Mediterranean cone snail, <i>Lautoconus ventricosus</i>. <a href="https://www.sciencedirect.com/science/article/pii/S1055790323001380">Molecular Phylogenetics and Evolution</a>.</li>
</ul>

---
