import time
import os
from dist_calc import distance

CURR_PATH = os.path.dirname(__file__)
os.chdir(CURR_PATH)

class StringUtils:
    def getOverlappingSubstrings(string: str, size: int) -> list[str]:
        substrings = [string[i:i+size] for i in range(0, len(string)-size+1)]
        return substrings

class MotifFinder:
    def __init__(self, N, d):
        self.maximum_mismatch = d
        self.sequence_count = N

        
class SequenceProcessor:
    class Lmer:
        def __init__(self, sequence, parent_idx, start_pos):
                self.sequence: str = sequence
                self.parent_idx = parent_idx
                self.start_pos = start_pos
                # lmer_map returns for each sequence number, the closest lmer from that sequence to this lmer if it is closer than 'd', otherwise it will return None
                self.instance_map = dict.fromkeys(list(range(Lmer.total_seq_count)), None)
                self.instance_map[seq_no] = (0, str.upper(lmer_seq))
                self.dist_buckets = [ [] for _ in range(Lmer.max_mismatch + 1) ]

        def __str__(self) -> str:
            return f'{self.parent_idx}:[{self.sequence}]'
        def __repr__(self) -> str:
            return self.__str__()
        
    def createLmer(sequence, parent_idx, start_pos)
    def readFromSitesFile(relative_file_path, n_seq = -1, seq_len = -1) -> str:
        unfiltered_strings = open(relative_file_path).read().lower().split('\n')
        sequences = [string[:seq_len] for string in unfiltered_strings if not string.startswith(">")]
        return sequences[:n_seq]
    
    def extractLmersFromSequence(sequence, lmer_size, sequence_no) -> list[Lmer]:
        lmer_list = []
        substrings = StringUtils.getOverlappingSubstrings(sequence, lmer_size)

        starting_pos = 0
        for substr in substrings:
            lmer_list.append( Lmer(substr, sequence_no, starting_pos) )
            starting_pos += 1

        return lmer_list
    
    def extractLmersFromSequenceList(string_list, size):
        seq_no = 0
        lmer_list = []

        for seq in string_list:
            lmer_list.append( SequenceProcessor.extractLmersFromSequence(seq, size, seq_no) )
            seq_no += 1

        return lmer_list
    
class BucketSearch:
    def fillBuckets(lmer_seq_list, d):

        lmer_count = len(lmer_seq_list)
        
        for i in range(lmer_count): # n
            for j in lmer_seq_list[i]: # t
                for k in range(i+1, lmer_count): # n
                    for l in lmer_seq_list[k]: # t -> n^2 t^2
                        dist = distance(j.lmer_seq, l.lmer_seq)
                        if dist <= d:
                            j.dist_buckets[dist].append(l)
                            l.dist_buckets[dist].append(j)

                            j_map = j.seq_map[k]
                            l_map = l.seq_map[i]

                            if j_map == None or (j_map != None and j_map[0] > dist):
                                j.seq_map[k] = (dist, l.lmer_seq)

                            if l_map == None or (l_map != None and l_map[0] > dist):
                                l.seq_map[i] = (dist, j.lmer_seq)
    
    def removeDegenerateLmers(lmer_seq_list):
        flat_list = [lmer for lmer_list in lmer_seq_list for lmer in lmer_list]
        i = 0
        while i < len(flat_list):
            lmer = flat_list[i]
            if None in lmer.seq_map.values():
                flat_list.pop(i)
            else: i += 1
        return flat_list
                
    def findMotif(candidate_lmers):
        least_distance = 999
        best_map = None

        for lmer in candidate_lmers:
            curr_total_dist = sum([t[0] for t in lmer.seq_map.values()])
            if curr_total_dist < least_distance:
                best_map = lmer.seq_map
                least_distance = curr_total_dist

        return best_map
            


def main():
    # sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0002.1.sites", n_seq=26)
    sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0014.2.sites", n_seq=30, seq_len=60)
    L = 10
    d = 6
    N = len(sequences)

    Lmer.total_seq_count = N 
    Lmer.max_mismatch = d
    
    start_time = time.perf_counter_ns()
    lmers = SequenceProcessor.extractLmersFromSequenceList(sequences, L)

    BucketSearch.fillBuckets(lmers, d)

    candidate_lmers = BucketSearch.removeDegenerateLmers(lmers)

    bestMotif = BucketSearch.findMotif(candidate_lmers)
    end_time = time.perf_counter_ns()

    # print(bestMotif[0].start_pos)
    # print(bestMotif)
    print(f"Time to find motif: {(end_time-start_time)/(10**6)} ms")

if __name__ == "__main__":
    main()