# DNA Toolset/Code testing file

from bio_seq import bio_seq

test_dna = bio_seq()
test_dna.generate_rnd_seq(length=300)

print(test_dna.get_seq_info())
#print(test_dna.get_bio_type())


print(test_dna.get_seq_info())
print(test_dna.gen_reading_frames())
print(test_dna.proteins_from_rf(["M","T","T","R","_"]))
print(test_dna.all_proteins())
