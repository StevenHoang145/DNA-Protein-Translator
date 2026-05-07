# DNA Bioinformatics Project

dna = input("Enter your DNA sequence: ").upper()

valid = ["A", "T", "C", "G"]

# Check if DNA is valid
good_dna = True

for base in dna:
    if base not in valid:
        good_dna = False

if good_dna == False:
    print("That is not a valid DNA sequence.")

else:

    print("\nDNA Sequence:")
    print(dna)

    print("\nLength:", len(dna))

    # Base counts
    a = dna.count("A")
    t = dna.count("T")
    c = dna.count("C")
    g = dna.count("G")

    print("\nBase Counts")
    print("A =", a)
    print("T =", t)
    print("C =", c)
    print("G =", g)

    # AT + GC counts
    at_count = a + t
    gc_count = g + c

    print("\nAT Count =", at_count)
    print("GC Count =", gc_count)

    # GC + AT percentages
    gc_percent = (gc_count / len(dna)) * 100
    at_percent = (at_count / len(dna)) * 100

    print("\nGC Content:", round(gc_percent, 2), "%")
    print("AT Content:", round(at_percent, 2), "%")

    # Melting temp
    melt_temp = 2 * at_count + 4 * gc_count

    print("Estimated Melting Temp:", melt_temp, "°C")

    # RNA sequence
    rna = dna.replace("T", "U")

    print("\nRNA:")
    print(rna)

    # Complement strand
    comp = ""

    for base in dna:

        if base == "A":
            comp += "T"

        elif base == "T":
            comp += "A"

        elif base == "C":
            comp += "G"

        elif base == "G":
            comp += "C"

    print("\nComplement:")
    print(comp)

    # Reverse stuff
    print("\nReverse DNA:")
    print(dna[::-1])

    print("\nReverse Complement:")
    print(comp[::-1])

    # Start codon
    if "ATG" in dna:
        print("\nStart codon found")

    # Stop codons
    stops = ["TAA", "TAG", "TGA"]

    for stop in stops:
        if stop in dna:
            print("Stop codon found:", stop)

    # Split into codons
    print("\nCodons:")

    codons = []

    for i in range(0, len(dna), 3):

        codon = dna[i:i+3]

        if len(codon) == 3:
            codons.append(codon)
            print(codon)

    # Codon usage
    print("\nCodon Counts")

    codon_count = {}

    for codon in codons:

        if codon in codon_count:
            codon_count[codon] += 1

        else:
            codon_count[codon] = 1

    for codon in codon_count:
        print(codon, ":", codon_count[codon])

    # Protein translation table
    table = {

        "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",
        "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
        "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
        "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",

        "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
        "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
        "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
        "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",

        "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
        "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
        "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
        "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",

        "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
        "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
        "TAC":"Y", "TAT":"Y", "TAA":"STOP", "TAG":"STOP",
        "TGC":"C", "TGT":"C", "TGA":"STOP", "TGG":"W"
    }

    print("\nProtein Translation")

    protein = ""

    for codon in codons:

        if codon in table:

            amino = table[codon]

            print(codon, "->", amino)

            if amino != "STOP":
                protein += amino + "-"

    print("\nProtein Sequence:")
    print(protein)

    # ORF finder
    print("\nLooking for ORFs...")

    found_orf = False

    for i in range(len(dna) - 2):

        if dna[i:i+3] == "ATG":

            for j in range(i + 3, len(dna) - 2, 3):

                current = dna[j:j+3]

                if current in stops:

                    print("Possible ORF:")
                    print(dna[i:j+3])

                    found_orf = True
                    break

    if found_orf == False:
        print("No ORFs found")

    # Restriction enzymes
    print("\nRestriction Enzymes")

    enzymes = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT"
    }

    found = False

    for enzyme in enzymes:

        site = enzymes[enzyme]

        if site in dna:

            print(
                enzyme,
                "site found at position",
                dna.find(site) + 1
            )

            found = True

    if found == False:
        print("No enzyme sites found")

    # Mutation comparison
    answer = input("\nCompare another sequence? ").lower()

    if answer == "yes":

        dna2 = input("Enter second DNA sequence: ").upper()

        print("\nMutation Results")

        smallest = min(len(dna), len(dna2))

        mutation = False

        for i in range(smallest):

            if dna[i] != dna2[i]:

                print(
                    "Position",
                    i + 1,
                    ":",
                    dna[i],
                    "->",
                    dna2[i]
                )

                mutation = True

        if len(dna2) > len(dna):
            print("Possible insertion mutation")

        elif len(dna2) < len(dna):
            print("Possible deletion mutation")

        if mutation == False:
            print("No mutations found")

    # Save file
    save = input("\nSave results to a text file? ").lower()

    if save == "yes":

        file = open("dna_results.txt", "w")

        file.write("DNA Analysis Results\n\n")

        file.write("DNA Sequence:\n")
        file.write(dna + "\n\n")

        file.write("Length: " + str(len(dna)) + "\n")

        file.write("A Count: " + str(a) + "\n")
        file.write("T Count: " + str(t) + "\n")
        file.write("C Count: " + str(c) + "\n")
        file.write("G Count: " + str(g) + "\n")

        file.write("AT Count: " + str(at_count) + "\n")
        file.write("GC Count: " + str(gc_count) + "\n")

        file.write("GC Content: " + str(round(gc_percent, 2)) + "%\n")
        file.write("AT Content: " + str(round(at_percent, 2)) + "%\n")

        file.write("RNA Sequence:\n")
        file.write(rna + "\n")

        file.write("Protein Sequence:\n")
        file.write(protein + "\n")

        file.close()

        print("Results saved to dna_results.txt")

    print("\nAnalysis complete")