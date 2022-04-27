 2\. Sequence_alignment 

Przyrównanie sekwencji
======================

    mkdir sequence_alignment
    cd sequence_alignment
    

Przyrównanie sekwencji to podstawowy problem genomiki, który jest częścią niemal każdych badań w tej dziedzinie. Przyrównując sekwencje określamy ich dopasowanie oraz podobieństwo. Idea przyrównania sekwencji biologicznych, oparta jest na założeniu, że sekwencje mające wspólne pochodzenie, wywodzące się z jednej sekwencji ancestralnej są do siebie bardziej podobne od sekwencji niespokrewnionych, i możliwe jest odtworzenie wzajemnej relacji odpowiadających sobie pozycji, nawet jeśli te uległy jakiejś zmianie (np. na drodze ewolucji). Przyrównywanie sekwencji pozwala nam na znalezienie podobieństw i różnic między genomami; określić które części ulegają zmianom, a które są raczej stabilne; zidentyfikowanie sekwencji powtórzonych; Generalnie wydzielić możemy dwa główne typy przyrównania sekwencji:

1.  Przyrównanie globalne - jest to przyrównanie co najmniej dwóch sekwencji od ich początków do końców, maksymalizujące pewną funkcję dopasowania. Algorytmem, który znajduje optymalne dopasowanie dwóch sekwencji jest algorytm Needlemana-Wunscha, który zaimplementowany jest np. w narzędziu [Emboss Needle](https://www.ebi.ac.uk/Tools/psa/emboss_needle/)
    
2.  Przyrównanie lokalne - przyrównanie znajdujące maksymalnie dopasowane fragmenty sekwencji. Tego typu przyrównanie możemy uzyskać wykorzystując algorytm Smitha-Watermana ([Emboss Water](https://www.ebi.ac.uk/Tools/psa/emboss_water/))
    

Ze względu na ilość i rozmiar analizowanych sekwencji, do przyrównań wykorzystuje się różne heurystyki, które chociaż nie zawsze generują optymalne przyrównania, są jednak zdecydowanie szybsze i mniej zasobochłonne. Powszechnie wykorzystywane metody przyrównania sekwencji oparte są o szybką identyfikację krótkich, dokładnych dopasowań (ziaren), rozszerzanie i weryfikację dopasowań. “Seed and extend” -> “Seed-chain-align”. Generalnie programy do przyrównywania sekwencji różnią się metodami wykorzystywanymi na każdym z tych etapów. Na każdym etapie przyrównania zoptymalizowanych może zostać wiele parametrów, w zależności od zastosowania oraz konkretnego przypadku.

W tym ćwiczeniu zapoznamy się z kilkoma klasycznymi sytuacjami, w których wykorzystujemy przyrównanie sekwencji, oraz poznamy przykładowe programy służące temu celu.

![figure 1](https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2Fs13059-021-02443-7/MediaObjects/13059_2021_2443_Fig1_HTML.png)  
[https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02443-7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02443-7)

### Przyrównanie całych genomów

Wygeneruj “genom” składający się z dwóch chromosomów o długości 100kbp każdy, zgodnie z HMM: A i B. Prawdopodobieństwo przejścia między stanami to: p(B|A) = 0.05, p(A|B) = 0.001. Prawdopodobieństwo emisji stanu A zgodne jest z rozkładem jendorodnym, prawdopodobieństwo emisji stanu B to:

|     | A   | C   | T   | G   |     |
| --- | --- | --- | --- | --- | --- |
| A   | 0.1 | 0.1 | 0.7 | 0.1 |     |
| T   | 0.7 | 0.1 | 0.1 | 0.1 |     |
| C   | 0.4 | 0.1 | 0.4 | 0.1 |     |
| G   | 0.4 | 0.1 | 0.4 | 0.1 |     |

Na potrzeby tego ćwiczenia możemy posłużyć się skryptem:

`python mo_zad1.py 100000 seq.fa`

Aby skomplikować nieco genom, wprowadzimy do niego trochę powtórzeń tandemowych oraz rozproszonych duplikacji. Wykorzystamy skrypt, który oprócz sekwencji DNA w formacie fasta przyjmuje ciąg operacji s - SNV; i - małe indele; l - duże indele; v - inwersje; t - duplikacje tandemowe; d - rozproszone duplikacje; n - translokacje niewzajemne (non-reciprocal); r - wzajemne translokacje.

    python variantize.py seq.fa tdsitdtdsi > seq1.fa
    

Wykorzystajmy profil kmerów do szybkiego porównania stopnia złożoności sekwencji:

    jellyfish count -m 21 -s 100M seq.fa -o 21mers.jf
    jellyfish histo 21mers.jf
    
    
    jellyfish count -m 21 -s 100M seq1.fa -o 21mers_2.jf
    jellyfish histo 21mers_2.jf
    

Wykorzystując program nucmer z pakietu [MUMmer](https://mummer4.github.io/manual/manual.html) dokonaj przyrównania dwóch genomów. Innymi często stosowanymi programami do przyrównania całych genomów są np. [BLAT](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC187518/), [LAST](https://gitlab.com/mcfrith/last), czy nowszy [GSAlign](https://github.com/hsinnan75/GSAlign).

    conda install mummer -c bioconda
    nucmer -h 
    

Przyrównanie sekwencji:

    nucmer --maxmatch seq.fa seq1.fa -p nucmer_out
    less nucmer_out.delta
    

Pzedstawienie w czytelnej formie:

    show-coords -Trlc nucmer_out.delta > nucmer_out.delta.coords
    

Wizualizacja za pomocą dotplotu:

    mummerplot  nucmer_out.delta
    mummerplot -x "[1000,10000]" -y ["1000,10000]" nucmer_out.delta
    

W pliku delta przechowywane są również informacje dotyczące przyrównań na poziomie nukleotydowym. Dzięki temu możemy zidentyfikować małe warianty, którymi różnią się przyrównywane genomy:

    show-snps 
    
    show-snps -T nucmer_out.delta|wc -l
    show-snps -TI nucmer_out.delta|wc -l
    show-snps -TCI nucmer_out.delta|wc -l
    

Aby ulepszyć identyfikację SNP możemy spróbować zidentyfikować najlepiej dopasowane powtórzenie:

    delta-filter -rq nucmer_out.delta > nucmer_out.filtered
    mummerplot nucmer_out.filtered
    show-snps -TCI nucmer_out.filtered|wc -l
    

Przyrównanie całego genomu możemy wykorzystać do identyfikacji sekwencji powtarzających się w obrębie jednego genomu:

    nucmer --maxmatch --nosimplify seq.fa seq.fa -p nucmer_out
    mummerplot  nucmer_out.delta
    
    nucmer --maxmatch --nosimplify seq1.fa seq1.fa -p nucmer_out
    mummerplot  nucmer_out.delta
    show-coords -rT nucmer_out.delta|less
    show-coords -rT nucmer_out.delta|awk '$1!=$3'
    

**Zadanie dodatkowe** \- Czy istnieje związek między kompresowalnością (LZ / gzip) a obecnością sekwencji powtórzonych?

#### Przyrównanie dwóch genomów przy pomocy [minimap2](https://github.com/lh3/minimap2) i wizualizacja za pomocą [dgenies](http://dgenies.toulouse.inra.fr/)

    conda install -c bioconda dgenies
    

Najlepiej odpalić w odrębnym terminalu:

    dgenies run  
    
    minimap2 -x asm5 seq.fa seq1.fa > aln.paf
    
    minimap2 -x asm5 -r100,1000 seq.fa seq1.fa > aln.paf
    
    minimap2 -x asm5 -r100,1000 seq1.fa seq1.fa  > aln.paf		
    
    minimap2 -x asm5 -r100,1000 -P seq1.fa seq1.fa  > aln.paf		
    
    minimap2 -x asm5 -r100,1000 -P seq.fa seq.fa  > aln.paf	
    

Szybki dotplot: [Emboss Dotmatcher](https://www.bioinformatics.nl/cgi-bin/emboss/dotmatcher)

**Zadanie**

Ze strony NCBI pobierz dwa genomy należące do tego samego gatunku (np. dwa szczepy bakterii) oraz do tego samego rodzaju (np. dwa gatunki bakterii). Przyrównaj genomy do siebie, zwizualizuj przyrównanie w formie dot-plotu. Zidentyfikuj i zintepretruj różnice między genomami. Czy pewne części genomu różnią się bardziej od pozostałych? Z czego może to wynikać? Przedstaw gęstość małych wariantów na genomie w formie wykresu liniowego.

#### Wyszukiwanie krótkich fragmentów homologicznych przy pomocy [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

    conda install blast -c bioconda
    
    jellyfish dump 21mers_2.jf|grep ^\>2 -A1|head -2 > query.fa
    
    
    blastn -subject seq1.fa -query query.fa
    

Domyślny program to megablast, z ziarnem o długości 28.  
[Parametry presetów blasta.](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/)

    blastn -subject seq1.fa -query query.fa -word_size 5
    
    blastn -subject seq1.fa -query query.fa -task blastn-short
    blastn -subject seq1.fa -query query.fa -task blastn-short -outfmt 6
    blastn -subject seq1.fa -query query.fa -task blastn-short -outfmt 7
    blastn -subject seq1.fa -query query.fa -task blastn-short -outfmt '6 qseqid qlen sseqid sstart send qstart qend evalue bitscore length pident mismatch gapopen sstrand sseq qseq'
    blastn -subject seq1.fa -query query.fa -task blastn-short -outfmt '6 qseqid qlen sseqid sstart send qstart qend evalue bitscore length pident mismatch gapopen sstrand sseq qseq'|awk '$8 < 1e-5'
    

Jeżeli nasze wyszukiwanie nie jest jednorazowe warto zbudować bazę danych sekwencji, korzystając z [makeblastdb](https://www.ncbi.nlm.nih.gov/books/NBK569841/)

### Mapowanie odczytów

Podstawą eksploracji genomów jest analiza ich sekwencji. Istnieje wiele sposobów sekwencjonowania genomów, jednak żaden z nich nie pozwala na odczytanie całej sekwencji. Zamiast tego otrzymywane są fragmenty o różnej długości i dokładności - w zależności od technologii. Kluczowym elementem wielu analiz genomowych jest zrekonstruowanie genomu, czyli odzyskanie informacji odnośnie wzajemnej relacji odczytów. Do tego celu wykorzystywane są dwa podejścia - składanie genomu (assembly) oraz - jeżeli istnieje złożenie, które możemy wykorzystać jako referencję - mapowanie odczytów do referencji. W tym ćwiczeniu dokonamy symulacji odczytów z wiodących technologii sekwencjonowania i zapoznamy się z metodyką mapowania odczytów do genomu na przykładzie kilku popularnych narzędzi. Dodatkowo zapoznamy się z metodami wizualizacji danych w przeglądarce genomowej IGV.

#### Illumina

[https://www.illumina.com/content/dam/illumina-marketing/documents/products/illumina\_sequencing\_introduction.pdf](https://www.illumina.com/content/dam/illumina-marketing/documents/products/illumina_sequencing_introduction.pdf)  
Obecnie największą zaletą Illuminy jest niska cena.

##### Symulacja odczytów

Zadanie 1. Na podstawie wygenerowanego genomu, korzystając z programu [ART](https://doi.org/10.1093/bioinformatics/btr708), dokonaj symulacji odczytów Illumina o następujących parametrach: platforma HiSeqX, biblioteka PCR free, odczyty sparowane, długość odczytu 150 bp, średnia długość fragmentu 350 bp, odchylenie standardowe średniej długości fragmentu 10 bp, trzydziestokrotne pokrycie. Następnie, po kontroli jakości i ewentualnym przygotowaniu, zmapuj odczyty do wyjściowego genomu przy pomocy bwa mem (domyślne parametry), posortuj oraz zindeksuj korzystając z [samtools](http://www.htslib.org/).

    conda install art -c bioconda
    
    art_illumina
    art_illumina -ss HSXt -sam -i seq1.fa -p -l 150 -f 30 -m 350 -s 10 -o illumina
    

##### Kontrola jakości

    conda install -c bioconda fastqc
    fastqc illumina*fq
    

##### Usuwanie adapterów

![enter image description here](https://supportassets.illumina.com/content/dam/illumina-support/images/bulletins/PEcell1.png)  
_[https://support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html](https://support.illumina.com/bulletins/2016/04/adapter-trimming-why-are-adapter-sequences-trimmed-from-only-the--ends-of-reads.html)_

W tym przypadku nie ma konieczności usuwania adapterów. Jeżeli jednak w naszych odczytach są pozostałości adapterów, możemy je usunąć przy pomocy np. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) czy [Cutadapt](https://cutadapt.readthedocs.io/en/stable/).

##### Mapowanie

    conda install -c bioconda bwa samtools
    

Tworzenie indeksu:

    bwa index seq1.fa
    

Mapowanie i sortowanie odczytów:

    bwa mem -t 10 seq1.fa illumina1.fq illumina2.fq|samtools sort > illumina_bwa.bam && samtools index illumina_bwa.bam
    

[https://samtools.github.io/hts-specs/SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf)

W drugiej kolumnie rekordu zakodowane są bitowe [Flagi](https://broadinstitute.github.io/picard/explain-flags.html), opisujące pewne właściwości zmapowanych odczytów.

##### Ocena jakości mapowania

Aby ocenić jakość mapowania, możemy wykorzystać programy w pakiecie samtools. Samtools flagstat podsumowuje statystyki flag:

    samtools flagstat
    samtools stats
    samtools stats -r seq1.fa illumina_bwa.bam|grep ^\#
    samtools stats -r seq1.fa illumina_bwa.bam|grep ^SN
    

Bardzo dobrej jakości raporty dotyczące mapowania generuje program [Qualimap](http://qualimap.conesalab.org/):

    conda install -c bioconda qualimap
    qualimap
    qualimap bamqc -bam illumina_bwa.bam -outdir qualimap_results
    

##### Filtrowanie odczytów

Samtools umożliwia nam odfiltrowanie odczytów  
samtools view illumina_bwa.bam “0:1-1000”|wc -l  
samtools view illumina_bwa.bam -f 2|wc -l  
samtools view illumina_bwa.bam -q 20|wc -l

Pytanie: W jaki sposób odfiltrować niezmapowane odczyty? Jakie może być znaczenie niezmapowanych odczytów?

##### Wizualizacja

Do wizualizacji zmapowanych odczytów wykorzystamy program [IGV](https://software.broadinstitute.org/software/igv/):

    wget https://data.broadinstitute.org/igv/projects/downloads/2.12/IGV_Linux_2.12.3_WithJava.zip 
    unzip IGV_Linux_2.12.3_WithJava.zip 
    IGV_Linux_2.12.3/igv.sh
    

### Nanopore

Odczyty Oxford Nanopore są obarczone dużo większym błędem, są jednak znacząco dłuższe. Obecny rekord to ponad 4 000 000 bp. Dodatkowo, zczytywany sygnał pochodzi z natywnych cząsteczek DNA, dzięki czemu możliwe jest wykrywanie modyfikacj DNA (np. metylacji)  
![https://www.nature.com/articles/s41587-021-01108-x/figures/1](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41587-021-01108-x/MediaObjects/41587_2021_1108_Fig1_HTML.png?as=webp)  
_[https://www.nature.com/articles/s41587-021-01108-x/figures/1](https://www.nature.com/articles/s41587-021-01108-x/figures/1)_

##### Symulacja odczytów

Do symulacji odczytów nanopore wykorzystamy narzędzie [Badread](https://github.com/rrwick/Badread):  
badread simulate --reference seq1.fa --quantity 20x > nanopore.fq

##### Mapowanie

Jednym z najbardziej popularnych narzędzi do mapowania danych z sekwencjonowania trzeciej generacji jest Minimap2, który poznaliśmy już w kontekście przyrównania całych genomów i wykorzystamy go równieżteraz . Innymi często wykorzystywanymi programami są: [NGMLR](https://github.com/philres/ngmlr), [LRA](https://github.com/ChaissonLab/LRA), [GraphMap](https://github.com/isovic/graphmap), [BLASR](https://github.com/mchaisso/blasr), [LAST](https://github.com/mcfrith/last-genome-alignments), [Mashmap](https://github.com/marbl/MashMap). Różnią się one czułością, specyficznością i szybkością. Istnieją również programy łączące zalety kilku maperów. Przykładowo [VULCAN](https://gitlab.com/treangenlab/vulcan) mapuje odczyty szybkim Minimap2, identyfikuje potencjalnie problematyczne regiony i remapuje je dokładniejszym NGMLR.

minimap2 -t 10 -ax map-ont seq1.fa nanopore.fq|samtools sort > nanopore.bam && samtools index nanopore.bam

Jeżeli to konieczne, adaptery możemy usunąć korzystając np z [Porechop](https://github.com/rrwick/Porechop).

### Pacbio HiFi

Odczyty HiFi stanowią kompromis między długością a jakością.  
![https://www.pacb.com/technology/hifi-sequencing/](https://www.pacb.com/wp-content/uploads/HiFi-reads-img.svg)  
![enter image description here](https://www.pacb.com/wp-content/uploads/img_hifi_reads.svg)  
_[https://www.pacb.com/technology/hifi-sequencing/](https://www.pacb.com/technology/hifi-sequencing/)_

##### Symulacja odczytów

badread simulate --reference seq1.fa --quantity 20x --error\_model pacbio2016 --qscore\_model ideal --identity 99.5,100,0.5 --start\_adapter\_seq “” --end\_adapter\_seq “” --length 10000,3000 > pacbio_hifi.fq

##### Mapowanie

minimap2 -t 10 -ax map-hifi seq1.fa pacbio\_hifi.fq|samtools sort > pacbio\_hifi.bam && samtools index pacbio_hifi.bam

### Przykład genomu diploidalnego

W praktyce rzadko zachodzi potrzeba mapowania odczytów do genomu, z którego pochodzą (chociaż jest to wykorzystywane podczas walidacji złożenia genomu). Generalnie odczyty mapujemy w celu zidentyfikowania różnic między genomami. Co więcej, w przypadku genomów eukariotycznych z reguły mamy do czynienia z wyższą ploidalnością, więc otrzymane odczyty pochodzą z różnych wersji chromosomów homologicznych. Wygenerujmy więc prosty przykład genomu diploidalnego, w którym spodziewać będziemy się zmienności w układzie homo- i heterozygotycznym.

python [variantize.py](http://variantize.py) seq1.fa si >seq2.fa  
python [variantize.py](http://variantize.py) seq2.fa si > seq3.fa  
cat seq2.fa seq3.fa > seq_diploid.fa

Zadanie 1. Na podstawie wygenerowanego genomu wygeneruj odczyty Illumina PE, Nanopore i HiFi, zgodnie z poprzednimi przykładami, oraz zmapuj je do genomu.

    art_illumina -ss HSXt -sam -i seq_diploid.fa -p -l 150 -f 15 -m 350 -s 10 -o illumina_diploid
    
    badread simulate --reference seq_diploid.fa --quantity 10x > nanopore_diploid.fq
    
    badread simulate --reference seq_diploid.fa --quantity 10x --error_model pacbio2016 --qscore_model ideal --identity 99.5,100,0.5 --start_adapter_seq "" --end_adapter_seq "" --length 10000,3000 > pacbio_hifi_diploid.fq
    
    bwa mem -t 10 seq1.fa illumina_diploid1.fq illumina_diploid2.fq|samtools sort > illumina_diploid_bwa.bam && samtools index illumina_diploid_bwa.bam
    
    minimap2 -t 10 -ax map-ont seq1.fa nanopore34.fq|samtools sort > nanopore34.bam && samtools index nanopore34.bam
    
    minimap2 -t 10 -ax map-hifi seq1.fa pacbio_hifi34.fq|samtools sort > pacbio_hifi34.bam && samtools index pacbio_hifi34.bam
    

Zadanie 2. Odfiltruj odczyty niezmapowane. Ile odczytów mapuje unikalnie?

Zadanie 3. Oblicz i porównaj statystyki flag (samtools flagstat) oraz qualimap dla wyników mapowań. Obejrzyj wyniki w IGV.

Zadanie 4. Wykorzystując samtools depth stworz wykres głębokości pokrycia odczytów w oknie 150, krok 10. Jaki procent genomu pokryty jest przynajmniej 5x? Czy są fragmenty do którego nie mapują żadne odczyty? Czy profil ulegnie zmianie, kiedy odfiltrujemy odczyty mapujące nieunikalnie?

Zadanie 5. Powtórz mapowanie krótkich odczytów, modyfikując długość ziarna do 2 i 150 bp. W jaki sposób zmiana wpływa na czas i wynik mapowania (ilość zmapowanych odczytów).

Zadanie Dodatkowe. Porównaj dwa dowolne mappery krótkich odczytów (np bwa-mem i [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml). Wykorzystując wiedzę na temat pochodzenia odczytów, oblicz ich precyzję (PPV, precision) p=TPTP+FPp = \\frac{TP}{TP+FP}p=TP+FPTP​ i czułość (True Positive Rate, sensitivity, recall) r=TPTP+FNr = \\frac{TP}{TP+FN}r=TP+FNTP​i wartość F (F-score) F1=2prp+rF1 = \\frac{2pr}{p+r}F1=p+r2pr​ Do wykonania zadania możesz wykorzystać program [Teaser](https://github.com/Cibiv/Teaser), lub napisać własną analizę.
