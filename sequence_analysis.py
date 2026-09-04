from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

# 1. Load TP53 CDS

sequence = """ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGA"""

print("TP53 CDS loaded successfully!")
print("Sequence length:", len(sequence))

print("A:", sequence.count("A"))
print("T:", sequence.count("T"))
print("C:", sequence.count("C"))
print("G:", sequence.count("G"))


# 2. Basic DNA analysis
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

gc = ((sequence.count("G") + sequence.count("C")) / len(sequence)) * 100

print("GC Content:", round(gc, 2), "%")


# 3. Codon analysis
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

first_codon = sequence[0:3]

print("First codon:", first_codon)
print("CDS divisible by 3:", len(sequence) % 3 == 0)
print("Number of codons:", len(sequence) // 3)


# 4. DNA → Protein translation
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

dna = Seq(sequence)
protein = dna.translate()

protein_clean = protein[0:-1]

print("Protein:", protein)
print("Protein length:", len(protein_clean), "aa")

print("Stop codon:", sequence[-3:])
print("Stop codon symbol:", protein[-1])


# 5. Reference protein validation
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

reference_protein = """MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"""

print("Reference protein length:", len(reference_protein), "aa")
print("Sequences are identical:", protein_clean == reference_protein)


# 6. Amino acid composition
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

amino_acids = [ "A", "R", "N", "D", "C","Q", "E", "G", "H", "I","L", "K", "M", "F", "P","S", "T", "W", "Y", "V" ]

print("\nAmino acid composition:")

for amino_acid in amino_acids:
    count = protein_clean.count(amino_acid)
    percentage = (count / len(protein_clean)) * 100

    print( amino_acid, "=", count, "(", round(percentage, 2), "%)" )

max_amino_acid = max( amino_acids, key=lambda aa: protein_clean.count(aa) )

min_amino_acid = min( amino_acids, key=lambda aa: protein_clean.count(aa) )

print( "Most frequent amino acid:", max_amino_acid, "=", protein_clean.count(max_amino_acid) )

print( "Least frequent amino acid:", min_amino_acid, "=", protein_clean.count(min_amino_acid) )
 

# 7. Protein molecular weight
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

mw = molecular_weight( protein_clean, seq_type="protein" )

print("Molecular weight:", round(mw, 2), "Da")
 

# 8. Example mutation analysis
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

original_codon = "GAG"
mutated_codon = "GTG"

print("\nMutation analysis:")
print("Original codon:", original_codon)
print("Mutated codon:", mutated_codon)

original_aa = Seq(original_codon).translate()
mutated_aa = Seq(mutated_codon).translate()

print("Original amino acid:", original_aa)
print("Mutated amino acid:", mutated_aa)


# 9. Final summary
print("//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////")

print("\n--- TP53 Analysis Summary ---")

print("CDS length:", len(sequence), "nt")
print("Number of codons:", len(sequence) // 3)
print("Protein length:", len(protein_clean), "aa")
print("GC Content:", round(gc, 2), "%")
print("Protein molecular weight:", round(mw, 2), "Da")

print( "Sequences are identical:", protein_clean == reference_protein )

print( "Most frequent amino acid:", max_amino_acid, "=", protein_clean.count(max_amino_acid) )

print( "Least frequent amino acid:", min_amino_acid, "=", protein_clean.count(min_amino_acid) ) 
