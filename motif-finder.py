import time
import os
from dist_calc import distance
from collections import defaultdict
from collections import Counter
from math import log2

CURR_PATH = os.path.dirname(__file__)
os.chdir(CURR_PATH)

class StringUtils:
    def getOverlappingSubstrings(string: str, size: int) -> list[str]:
        substrings = [string[i:i+size] for i in range(0, len(string)-size+1)]
        return substrings

    def displayMotifs(sequences, instances, motif_prototype):
        result = ""
        first_instance = {i: instances[i][0][1] for i in instances}

        result += f'{motif_prototype.seq_idx}\t'
        result += sequences[motif_prototype.seq_idx][:motif_prototype.start_pos]
        result += '\033[92m' + sequences[motif_prototype.seq_idx][motif_prototype.start_pos:motif_prototype.start_pos+10] + '\033[0m'
        result += sequences[motif_prototype.seq_idx][motif_prototype.start_pos+10:] + '\n'
    
        for seq_idx, start_pos in first_instance.items():
            result += f'{seq_idx}\t'
            result += sequences[seq_idx][:start_pos]
            result += '\033[91m' + sequences[seq_idx][start_pos:start_pos+10] + '\033[0m'
            result += sequences[seq_idx][start_pos+10:] + '\n'
        
        print(result)

    def IC(sequences: list[str]):
        information_content = 0
        string_length = len(sequences[0])
        seq_count = len(sequences)
        concat_string = ''.join(sequences)
        counts = Counter(concat_string)
        A_total = counts['A']
        C_total = counts['C']
        G_total = counts['G']
        T_total = counts['T']

        for p in range(string_length):
            A_odds = sum([c=='A' for s in sequences for c in s[p]])/seq_count
            C_odds = sum([c=='C' for s in sequences for c in s[p]])/seq_count
            G_odds = sum([c=='G' for s in sequences for c in s[p]])/seq_count
            T_odds = sum([c=='T' for s in sequences for c in s[p]])/seq_count
            information_content += log2(A_odds/A_total) + log2(C_odds/C_total) + log2(G_odds/G_total) + log2(T_odds/T_total)

        return information_content


        
class SequenceProcessor:
    def __init__(self, N, w):
        self.maximum_wiggle_room = w
        self.sequence_count = N
    class Lmer:
        def __init__(self, lsequence, seq_idx, start_pos, w, N):
                self.lsequence: str = lsequence
                self.seq_idx = seq_idx
                self.start_pos = start_pos

                self.instances = {k: SequenceProcessor.Instance(w) for k in range(N)}
                del self.instances[self.seq_idx]

                self.instance_mask = [False for _ in range(N)]
                self.instance_mask[self.seq_idx] = True
        def __str__(self) -> str:
            return f'{self.seq_idx}:[{self.lsequence}]'
        def __repr__(self) -> str:
            return self.__str__()
    
    class Instance:
        def __init__(self, w) -> None:
            self.w = w
            self.least_distance = 999
            self.instance_dict = defaultdict(list)
        def appendIfNear(self, dist, lmer):
            if dist <= self.least_distance + self.w: # Append Lmer to instances
                self.instance_dict[dist].append(lmer)
                if dist < self.least_distance: # If this new Lmer is also the new closest, remove some of the older ones that are too far.
                    self.least_distance = dist
                    keys = list(self.instance_dict.keys())
                    for key in keys:
                        if key - dist > self.w:
                            del self.instance_dict[key]
        def __str__(self) -> str:
            return f'{self.instance_dict}'
        def __repr__(self) -> str:
            return self.__str__()
        
    def createLmer(self, lsequence, seq_idx, start_pos):
        return self.Lmer(lsequence, seq_idx, start_pos, self.maximum_wiggle_room, self.sequence_count)
    
    def readFromSitesFile(relative_file_path, n_seq = -1, seq_len = -1) -> str:
        unfiltered_strings = open(relative_file_path).read().upper().split('\n')
        sequences = [string[:seq_len] for string in unfiltered_strings if not string.startswith(">")]
        sequences = [sequences[i] for i in range(0, len(sequences), 2)]
        return sequences[:n_seq]
    
    def extractLmersFromSequence(self, sequence, seq_idx, lmer_size) -> list[Lmer]:
        lmer_list = []
        lmer_sequences = StringUtils.getOverlappingSubstrings(sequence, lmer_size)

        starting_pos = 0
        for lseq in lmer_sequences:
            lmer_list.append( self.createLmer(lseq, seq_idx, starting_pos) )
            starting_pos += 1

        return lmer_list
    
    def extractLmersFromSequenceList(self, sequence_list, size):
        seq_idx = 0
        lmer_list = []

        for seq in sequence_list:
            lmer_list.append( self.extractLmersFromSequence(seq, seq_idx, size) )
            seq_idx += 1

        return lmer_list
    
class MotifFinder:
    def findInstances(lmer_seq_list, N, d):
        for first_seq in range(N):
            for first_lmer in lmer_seq_list[first_seq]:
                for second_seq in range(first_seq+1, N):
                    for second_lmer in lmer_seq_list[second_seq]:
                        dist = distance(first_lmer.lsequence, second_lmer.lsequence)
                        if dist <= d:
                            first_lmer.instance_mask[second_seq] = True
                            first_lmer.instances[second_seq].appendIfNear(dist, second_lmer)

                            second_lmer.instance_mask[first_seq] = True
                            second_lmer.instances[first_seq].appendIfNear(dist, first_lmer)

    
    def removeDegenerateLmers(lmer_list, seq_count):
        candidates = []
        for lmer in lmer_list:
            if sum(lmer.instance_mask) == seq_count:
                candidates.append(lmer)
        return candidates
    
    def findBestMotif(candidate_prototypes, seq_count):
        best_motif_score = 0

        for prototype in candidate_prototypes:
            instance_list = [prototype.lsequence]
            other_sequences = list(range(seq_count)).remove(prototype.seq_idx)
            for i in other_sequences:
                
    

def main():
    # sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0002.1.sites", n_seq=26)
    sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0014.2.sites", n_seq=30, seq_len=60)
    L = 10
    d = 5
    w = 1
    N = len(sequences)
    sq = SequenceProcessor(N, w)
    
    start_time = time.perf_counter_ns()
    ########
    lmer_seq_lst = sq.extractLmersFromSequenceList(sequences, L)
    lmer_flatlist = [lmer for lmer_list in lmer_seq_lst for lmer in lmer_list]

    MotifFinder.findInstances(lmer_seq_lst, N, d)
    candidate_lmers = MotifFinder.removeDegenerateLmers(lmer_flatlist, N)
    # print(len(candidate_lmers))
    # motif = MotifFinder.findMotif(candidate_lmers)
    ########
    end_time = time.perf_counter_ns()

    # StringUtils.displayMotifs(sequences, motif[0], motif[1])
    print(f"Time to find motif: {(end_time-start_time)/(10**6)} ms")

if __name__ == "__main__":
    main()