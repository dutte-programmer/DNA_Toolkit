from bio_structs import *
import random
import collections

class bio_seq:
    ''' DNA sequence class. Default value: ATCG, DNA, No label '''

    def __init__(self, seq="ATCG", seq_type="DNA", label="No_Label"):
        '''Sequence initialization, valdiations'''
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, "Provided data does not seem to be correct {}".format(self.seq_type)

    def __validate(self):
        '''Validates if Sequence only contains ACGT '''
        return set(Nucleotides).issuperset(self.seq)

    def get_seq_info(self):
        '''Returns 4 strings. Full sequence information'''
        return "[Label]: {}\n[Sequence]: {}\n[Biotype]: {}\n[Length]: {}".format(self.label, self.seq, self.seq_type, len(self.seq))

    def get_bio_type(self):
        '''Returnst biotype'''
        return self.seq_type

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        '''Generate a random sequence, provided the length'''
        seq = "".join([random.choice(Nucleotides) for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")


    def count_nuc_frequencies(self):
        frequencies = {"A": 0, "C": 0, "G": 0, "T": 0}
        for nuc in self.seq:
            frequencies[nuc] += 1
        return frequencies
        # alternative optimization: return dict(collections.Counter(seq))

    def transcription(self):
        temporary = self.seq.lower()
        tmp1 = temporary.replace("a", "U")
        tmp2 = tmp1.replace("t", "A")
        tmp3 = tmp2.replace("c", "G")
        tmp4 = tmp3.replace("g", "C")
        return tmp4

    def complement(self):
        complement1 = "".join([DNA_reversecomplement[nuc] for nuc in self.seq])
        return complement1

    def reverse_complement(self):
        reverse = self.seq[::-1]
        complement1 = "".join([DNA_reversecomplement[nuc] for nuc in reverse])
        return complement1

    def GC_content(self):
        '''DNA/RNA GC Count'''
        return round(self.seq.count("C") + self.seq.count("G") / len(self.seq) * 100)

    def GC_content_subsec(self, k=20):
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round(subseq.count("C") + subseq.count("G") / len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        '''Translates a sequence into a single letter Aminoacid Code'''
        return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        ''' Provides the frequency of each codon encoding a given aminoacid in a DNA sequence'''
        temporary = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                temporary.append(self.seq[i:i + 3])

        freqDict = dict(collections.Counter(temporary))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        '''Generate the six reading frames of a DNA sequence, including the reverse complement'''
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self, aa_seq):
        proteins = []
        current_prot = []
        for aa in aa_seq:
            if aa == "_":
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return [x + "*" for x in proteins]

    def all_proteins(self, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res