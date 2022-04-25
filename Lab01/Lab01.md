  1\. Genome\_characteristics 

Genom - sekwencja nukleotydów
=============================

### Wprowadzenie

W genomach zgromadzona jest większość informacji niezbędnej do funkcjonowania komórek i zestaw instrukcji niezbędnych do powielenia samych siebie. “Instrukcje” zawarte na polimerze kwasu deoksyrybonukleinowego warunkują strukturę i organizację przestrzenną DNA; określają sekwecje, a zatem również struktury i funkcje RNA i białek oraz stanowią podstawę regulalacji powstawania tych cząsteczek. Na elementarnym poziomie ogrom różnorodnej informacji kodowanej przez genomy określony jest przez sekwencje złożone z czteroliterowego alfabetu {A,C,T,G}. Chociaż obecna analiza genomów opiera się przede wszystkim o dotychczas zgromadzoną wiedzę odnośnie znaczenia poszczególnych “słów”, informacja zawarta w surowej sekwencji - jej skład i ogólny stopień i charakter uporządkowania pozwala nam zdobyć pewien szeroki ogląd na genom, dostrzec jego niejednorodność i zauważyć sygnały mogące mieć znaczenie biologiczne.

### Cel zajęć

W tym ćwiczeniu scharakteryzujemy genomy korzystając z naiwnych, nie odnoszących się do konkretnej wiedzy biologicznej miar, zastosujemy statystyczne metody pozwalające na ocenę przypadkowości obserwowanych sygnałów, oraz zwizualizujemy uzyskane wyniki.

### Środowisko

W trakcie ćwiczeń korzystać będziemy z systemu zarzadzania pakietami [Conda](https://conda.io). Sugeruję stworzyć jedno środowisko, do którego będziemy wracać w kolejnych ćwiczeniach. Niektóre programy będą wymagały stworzenia osobnego środowiska, ze wzgldu na niekompatybilność wersji wymaganych pakietów.

Instalacja Condy:
```bash
    wget  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    
    chmod -x ./Miniconda3-latest-Linux-x86_64.sh
    
    ./Miniconda3-latest-Linux-x86_64.sh
    
    /home/student/miniconda3/bin/conda init
    
    source ~/.bashrc
    
    conda config --add channels bioconda 
    
    conda config --add channels conda-forge
    
    conda install mamba
 ```   

Stwórz środowisko w conda:
```bash
    conda create -n genomika22 -c bioconda
    conda activate genomika22
    conda install python=3
```
    

Utwórz katalog roboczy
```bash
    mkdir genome_statistics
    cd genome_statistics
```

### Probabilistyczne modele sekwencji genomowych

W tym ćwiczeniu, podobnie jak w części kolejnych zajęć opierać będziemy się na symulacjach sekwencji genomowych, o rosnącym stopniu złożoności. Wiele problemów napotykanych w trakcie analiz genomicznych, takich jak dzisiejsze statystyki opisowe, przyrównywanie sekwencji, składanie genomów czy identyfikacja wariantów zależna jest jedynie od samej sekwencji, nie od jej znaczenia biologicznego. Wykorzystanie symulacji pozwoli nam na odniesienie otrzymanych w wyników do znanego nam modelu wyjściowego, i tym samym ocenę przeprowadzonej analizy. Dzięki temu poznamy wzajemny wpływ cech sekwencji i parametrów poszczególnych analiz. Istotną kwestią jest również zużycie zasobów obliczeniowych. Rozmiary genomów bywają bardzo duże, a ich analizy kosztowne i czasochłonne. Symulacje pozwolą na przećwiczenie zróżnicowanych podejść w trakcie zajęć, a poznane techniki będą mogły zostać wykorzystane do analizy realnych danych poza nimi.

#### Modele wielomianowe

Najprostszy model sekwencji zakłada niezależność pozycji i równomierne rozmieszczenie znaków. Sekwencja powstaje w wyniku procesu stochastycznego w którym znaki są kolejno losowane z alfabetu {A,C,T,G}. Prawdopodobieństwo wylosowania poszczególnych znaków może, lecz nie musi być jednostajne. Modele wielomianowe często traktowane są jako model zerowy, odstępstwa od którego mogą wskazywać na pewne interesujące regiony genomu.

**Zadanie 1.** Wygeneruj losową sekwencję nukleotydów o długości 1000 bp przy założeniu jednorodnego (jednostajnego) rozkładu nukleotydów. Wynik zapisz w formacie fasta.
```python
    import random
    nucleotides = ['A', 'T', 'C', 'G']
    length = 1000
    
    def random_seq(nucleotides, n):
    	random_set = [random.choice(nucleotides) for i in range(n)]
    	return ''.join(random_set)
    
    seq = random_seq(nucleotides, length)
    print(seq)
```    

**Zadanie 2.** Wygeneruj losową sekwencję nukleotydów o długości 1000 bp zakładając niejednorodny, arbitralny rozkład nukleotydów p(A) = 0.1, p(C) = 0.3, p(T) = 0.2. Wynik zapisz w formacie fasta.
```python
    weights = [0.1, ...]
```    

albo:
```python
    nucleotides = {'A': 0.1, ...}
```    

na przykład:
```python
    random.choices(list, weights, k)
```    

#### Modele Markova jako modele generatywne.

Klasą modeli często wykorzystywanych do analizy danych genomowych są modele Markova. Modele Markova od modeli wielomianowych (Procesu Bernouliego) odróżnia brak niezależności - prawdopodobieństwo zaobserwowania symbolu S zależne jest w pewien sposób od poprzedzającej go sekwencji.

Najprostrzymi modelami Markova są łańcuchy Markova (Markov Chains), w przypadku których prawdopodobieństwo zaobserwowania symbolu S(n) na pozycji n zależy od sekwencji S(n-k:n-1) dla k-rzędowego łańcucha.

**Zadanie 3.** Wygeneruj losową sekwencję nukleotydów o długości 1000 bp zgodnie ze zilustrowanym na diagramie I-rzędowym procesem Markowa, zakładając jednorodny rozkład początkowy. Wynik zapisz w formacie fasta.  
![Image caption](https://i2.wp.com/1.bp.blogspot.com/-_v4xzlnHqWk/T4jy6aaBEhI/AAAAAAAAAo8/_L0D5t0VS7o/s400/first.png?resize=400%2C347)
```python
    transition_matrix = {
    'A': {'A':0.6, 'C':0.1, 'T':0.2, 'G':0.1},
    'C': {'A':0.1, 'C':0.5, 'T':0.1, 'G':0.3},
    'T': {'A':0.4, 'C':0.05, 'T':0.5, 'G':0.05},
    'G': {'A':0.05, 'C':0.2, 'T':0.05, 'G':0.7}
    }
```    

**Zadanie 4.** W jaki sposób należałoby zmodyfikować diagram, aby był równoważny z modelami z zadania 1 i 2?

Bardziej złożonymi modelami są Ukryte Modele Markova (HMM - Hidden Markov Models). W ich przypadku sekwencja definiowana jest przez zestaw stanów i prawdopodobieństw tranzycji między nimi, gdzie każdy stan cechuje inny rozkład prawdopodobieństwa emisji symboli. Modele te znalazły szerokie zastosowanie do modelowaniu zróżnicowanych problemów biologicznych - np. segmentacji sekwencji (np. egzon/intron), predykcji genów czy modelowania architektury sekwencji białkowych. Dobrym wprowadzeniem do HMM w analizie sekwencji biologicznych jest krótki artykuł:

> Eddy, S. What is a hidden Markov model?. _Nat Biotechnol_ **22,** 1315–1316 (2004). [https://doi.org/10.1038/nbt1004-1315](https://doi.org/10.1038/nbt1004-1315)

Z HMM spotkamy się jeszcze podczas zajęć dotyczących adnotacji sekwencji biologicznych, dzisiaj wykorzystamy je w roli modeli generatywnych.

**Zadanie 5.** Wygeneruj losową sekwencję nukleotydów o długości 1000 bp, składającą się z dwóch stanów ukrytych: A i B. Prawdopodobieństwo przejścia między stanami to: p(B|A) = 0.05, p(A|B) = 0.1. Prawdopodobieństwo emisji nukleotydów ze stanu A zgodne jest z rozkładem jednorodnym, prawdopodobieństwo emisji nukleotydów ze stanu B zgodne jest z rozkładem wielomianowym zdefiniowanym w Zadaniu 2. Załóż jednorodny rozkład początkowy. Wynik zapisz w formacie fasta.

### Charakterystyka sekwencji

W tej części scharakteryzujemy nasze"genomy" korzystając z kilku miar kompozycji oraz złożoności sekwencji. Zastosujemy je dla całego genomu, oraz w oknie. W ten sposób zlokalizujemy fragmenty które odbiegają od rozkładu losowego i w przypadku realnych danych mogą być istotne biologicznie.

**Zadanie 6.** Dla każdej z wygenerowanych sekwencji oblicz:

1.  Relatywny udział poszczególnych nukleotydów
2.  Zawartość GC (GC content)
3.  Asymetrię GC i AT (GC / AT skew) - stosunek G do C (A do T) na danej nici DNA
4.  Entropia Shannona (miara nieuporządkowania sekwencji, czyli zawartości informacji)  
    H(X)\=−∑i\=1nP(xi)logP(xi)H(X) = -\\sum^n\_{i=1}P(x\_i)logP(x\_i)H(X)\=−i\=1∑n​P(xi​)logP(xi​)
5.  Złożoność Lempel’a-Ziv’a - jest to liczba unikalnych podsłów (podsekwencji) spotykanych podczas czytania sekwencji od lewej do prawej. Może być wykrzystana jako miara powtarzalności sekwencji. Przykładowo: sekwencję ACTGTGATCCTGACTGA rozłożymy do A|C|T|G|TGA|TC|CTGA

Obliczenia dotyczące kompozycji sekwencji możemy wykonać na przykład za pomocą zbioru narzędzi [seqkit](https://bioinf.shenwei.me/seqkit/)
```bash
    conda install seqkit -c bioconda
```    

Udział A:
```bash
    seqkit fx2tab in.fa -n -B A
```    

Zawartość GC:
```bash
    seqkit fx2tab in.fa -n -B GC
    seqkit fx2tab in.fa -n -g
```    

Asymetria GC:
```bash
    seqkit fx2tab in.fa -n -g
```    

**Zadanie 7.** Powtórz obliczenia z zadania 6 stosując okno 20 bp i krok 5 bp. Zwizualizuj wyniki na wykresie liniowym. Opisz otrzymany wynik.
```bash
    conda install matplotlib pandas
```    

GC w seqkit:
```bash
    cat in.fa| \
    seqkit sliding -s 5 -W 20 | \
    seqkit fx2tab -n -g > gc_sliding
```    

### Profile k-merów oraz ocena ich losowości

Profil k-merów przedstawia zliczenia wszystkich podsłów o długości k danej sekwencji. Znając profil k-merów możemy zidentyfikować nad- i niedoreprezentowane sekwencje, z których obie klasy mogą być istotne biologicznie. Przykładowo słów determinujących miejsca wiązań czynników transkrypcyjnych możemy spodziewać się rzadko i tylko w określonych miejscach.

**Zadanie 8.** Dla sekwencji 1, 2 i 3 oblicz obserwowane częstości 2-merów oraz określ teoretyczne prawdopodobieństwo ich wystąpienia przy założeniu modelu wielomianowego. Oblicz iloraz szans i logarytm ilorazu szans. Opisz otrzymany wynik.

P(s)\=P(s1s2...sn)P(s) = P(s\_1s\_2...s\_n)P(s)\=P(s1​s2​...sn​)  
or\=N(xy)N(x)N(y)or = \\frac{N(xy)}{N(x)N(y)}or\=N(x)N(y)N(xy)​

Powszechnie wykorzystywanym programem do analizy k-merów jest [jellyfish](https://github.com/gmarcais/Jellyfish)
```bash
    conda install jellyfish -c bioconda
```    

2 mery:
```bash
    jellyfish count -m 2 -s 100M -t 10 -C in.fa -o twomers.jf
```    

Aby wyświetlić zawartość pliku:
```bash
    jellyfish dump twomers.jf
```    

Liczebność AA:
```bash
    jellyfish query twomers.jf AA
```    

Histogram (profil k-merów):
```bash
    jellyfish histo twomers.jf
```    

Statystyki:
```bash
    jellyfish stats twomers.jf 
```    

Jednym ze sposobów na sprawdzenie losowości zaobserowanego sygnału często stosowanym w analizie sekwencji są randomizacje i metody Monte Carlo. Najprostszym przykładem są permutacje, czyli losowanie bez zwracania. Z sekwencji o długości L losujemy, bez zwracania nową sekwencję o długości L, dla której obliczamy interesującą nas statystykę (np. częstość występowania określonego k-meru). Proces powtarzamy N razy. W ten sposób otrzymujemy rozkład prawdopodobieństwa krotności zaobserwowania danego k-meru w losowej sekwencji o takich samych parametrach. Inną często stosowaną metodą jest bootstrap, w przypadku której losujemy ze zwracaniem (innymi słowy - losujemy z rozkładu wielomianowego z empirycznie określonymi częstościami występowania nukleotydów). W ten sposób możemy modelować populację, z której pochodzi nasza próba.

**Zadanie 9.** Wykorzystując test permutacyjny (100 iteracji, jednonukleotydowe próbkowanie) oszacuj rozkład częstości najczęściej spotykanego dinukleotydu w sekwencji 3. Określ prawdopodobieństwo wystąpienia częstości o wartości równej zaobserwowanej lub wyższej. Zobrazuj wynik na wykresie.

### Kolista wizualizacja danych genomowych - biblioteka Circos

Często wykorzystywanym sposobem wizualizacji danych genomowych są wykresy oparte o okrąg, rozpropagowane w dużej mierze przez Martina Krzywinskiego, twórcę biblioteki [Circos](http://circos.ca/). Dużą zaletą takiej reprezentacji jest możliwośc zwizualizowania zależności (np. regionów syntenicznych, czy interakcji) między fragmentami genomów. Pakiet Circos jest niezwykle elastyczny i doskonale udokumentowany, [z wieloma przykładami i dokładnymi opisami funkcjonalności.](http://circos.ca/documentation/tutorials/quick_start/) W tym ćwiczeniu stworzymy prosty wykres cyrkularny jednego z genomów, na którym zwizualizujemy uzyskane dotychczas wyniki.

![http://circos.ca/intro/genomic_data/](http://circos.ca/intro/genomic_data/img/circos-conde-nast-large.png)  
_[http://circos.ca/intro/genomic\_data/](http://circos.ca/intro/genomic_data/)_

**Zadanie 10.** Wizualizacja danych genomowych z pakietem Circos
```bash
    conda install circos -c bioconda
    circos -conf config.conf
```    

Aby stworzyć wykres, musimy przygotować kilka plików konfiguracyjnych:

1.  [Kariotyp](http://circos.ca/documentation/tutorials/ideograms/karyotypes/) - plik, w którym określamy podstawowe informacje dotyczące segmentów (chromosomów), które chcemy zwizualizować. Format:

| **chr** | **-** | **ID (w danych)** | **etykieta** | **początek** | **koniec** | **kolor** |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| chr | - | 1 | chromosome1 | 0 | 580076 | purple |


2.  [Główny plik konfiguracyjny](http://circos.ca/documentation/tutorials/quick_start/hello_world/) - plik, w którym definiujemy elementy grafiki, oraz ich parametry. W głównym pliku konfiguracyjnym określamy lokalizację pliku z kariotypem:
```html
    karyotype = karyo.conf
```
4.  Elementami, które możemy uwzględnić na wykresie są [m.in](http://m.in).:

Ideogramy, czyli segmenty “chromosomów”:

```html
    <ideogram>
        <spacing> default = 0.005r </spacing>
      
        radius = 0.90r  
        thickness = 20p  
        fill = yes
      
        #stroke\_thickness = 1  
        #stroke\_color = black
    </ideogram>
```

Oznaczenia osi (ticks):

```html
    <tick> spacing = 10000 u  
        color = grey  
        size = 10p  
    </tick>
```

Wykresy prezentujące dane liczbowe:

```html
    <plot>
        type = histogram  
        file = gc\_content.histo  
        thickness = 0p
    </plot>
```

Linie łączące elementy wykresu:

```html
    <links>  
        <link> radius = 0.8r  
        bezier\_radius = 0r  
        bezier\_radius\_purity = 0.9  
        color = black  
        thickness = 2  
        file = max\_kmer.links  
        </link>  
    </links>
```

5.  Dodatkowe pliki konfiguracyjne - możemy definiować w nich te same elementy, co w głównym pliku konfiguracyjnym. Uwzględniamy je w głównym pliku linijką <<include name.conf>>.
    
6.  Pliki z danymi do wykresów w formacie:

| **chr** | **start** | **stop** | **wartość** | **kolor** |
| :---: | :---: | :---: | :---: | :---: |
| 0 | 0 | 100 | 1.629407336396254 | fill_color=red |
| 0 | 100 | 200 | 1.5523160283080437 | fill_color=red |

7.  Pliki z danymi do połączeń w formacie:

| **id** | **chrom** | **start** | **stop** | **kolor** |
| :---: | :---: | :---: | :---: | :---: |
| 01 | 0 | 169475 | 169522 | color=blue |
| 01 | 0 | 224535 | 224555 | color=blue |

Wykres tworzymy poleceniem:
  
```bash
    circos -conf circos_conf.conf
```    

### Zadanie domowe:

Z bazy danych NCBI ( [https://www.ncbi.nlm.nih.gov/datasets/genomes/](https://www.ncbi.nlm.nih.gov/datasets/genomes/)) pobierz wybrany przez siebie genom referencyjny, opublikowany nie wcześniej niż w 2018 roku. Scharakteryzuj go pod względem miar poznanych na zajęciach. Utwórz profil 21-merów. Jakie jest prawdopodobieństwo, że zaobserwowana liczebność najliczniejszego 21-meru jest przypadkowa? Zwizualizuj wyniki na kolistym wykresie.

* * *
