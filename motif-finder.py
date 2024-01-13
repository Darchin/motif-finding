import time
import os
from dist_calc import distance

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

        
class SequenceProcessor:
    def __init__(self, N, d):
        self.maximum_mismatch = d
        self.sequence_count = N
    class Lmer:
        def __init__(self, lsequence, seq_idx, start_pos, d, N):
                self.lsequence: str = lsequence
                self.seq_idx = seq_idx
                self.start_pos = start_pos
                self.instance_map = {k: [] for k in range(N)}
                self.closest_instance_dist = {k: 999 for k in range(N)}
                del self.instance_map[self.seq_idx]
                del self.closest_instance_dist[self.seq_idx]
                self.dist_buckets = [[] for _ in range(d + 1)]

        def __str__(self) -> str:
            return f'{self.seq_idx}:[{self.lsequence}]'
        def __repr__(self) -> str:
            return self.__str__()
        
    def createLmer(self, lsequence, seq_idx, start_pos):
        return self.Lmer(lsequence, seq_idx, start_pos, self.maximum_mismatch, self.sequence_count)
    
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
    def fillDistanceBuckets(lmer_seq_list, N, d, w):
        for i in range(N): # n
            for j in lmer_seq_list[i]: # t
                for k in range(i+1, N): # n
                    for l in lmer_seq_list[k]: # t -> n^2 t^2
                        dist = distance(j.lsequence, l.lsequence)
                        if dist <= d:
                            j.dist_buckets[dist].append(l)
                            l.dist_buckets[dist].append(j)
                            if dist < j.closest_instance_dist[k]: j.closest_instance_dist[k] = dist
                            if dist < l.closest_instance_dist[i]: l.closest_instance_dist[i] = dist
    def findInstances(lmer_list, d, w):
        for lmer in lmer_list: #nt
            for i in range(d+1): # every distance up to d mismatch -> d
                for lmers_i_dist_away in lmer.dist_buckets[i]: # lmers 'i' distance away from 'lmer'
                    for lmer_at_i_dist in lmers_i_dist_away: # each lmer
                        closest = lmer.closest_instance_dist[lmer_at_i_dist.seq_idx] # closest_instance_from_same_seq
                        if i <= closest + w:
                            lmer.instance_map[lmer_at_i_dist.seq_idx].append(lmer_at_i_dist)

    
    def removeDegenerateLmers(lmer_list):
        candidates = lmer_list
        i = 0
        while i < len(lmer_list):
            # removed = False
            if 999 in candidates[i].closest_instance_dist.values():
                candidates.pop(i)
            else:
                i += 1
            # for val in lmer_list[i].instance_map.values():
            #     if len(val) == 0:
            #         removed = True
            #         lmer_list.pop(i)
            # if not removed:
                # i += 1

        return candidates
                
    def findMotif(candidate_lmers):
        least_distance = 999
        best_instances: dict
        motif_prototype = 0
        for lmer in candidate_lmers:
            curr_total_dist = sum([t[0][0] for t in lmer.instance_map.values()])
            if curr_total_dist < least_distance:
                best_instances = lmer.instance_map
                motif_prototype = lmer
                least_distance = curr_total_dist

        return best_instances, motif_prototype
    
        

def main():
    # sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0002.1.sites", n_seq=26)
    sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0014.2.sites", n_seq=30, seq_len=60)
    L = 10
    d = 5
    w = 1
    N = len(sequences)

    sq = SequenceProcessor(N, d)
    
    start_time = time.perf_counter_ns()
    ########
    lmer_seq_lst = sq.extractLmersFromSequenceList(sequences, L)
    lmer_flatlist = [lmer for lmer_list in lmer_seq_lst for lmer in lmer_list]

    MotifFinder.fillDistanceBuckets(lmer_seq_lst, N, d)


    candidate_lmers = MotifFinder.removeDegenerateLmers(lmer_flatlist)

    motif = MotifFinder.findMotif(candidate_lmers)
    ########
    end_time = time.perf_counter_ns()

    StringUtils.displayMotifs(sequences, motif[0], motif[1])
    print(f"Time to find motif: {(end_time-start_time)/(10**6)} ms")

if __name__ == "__main__":
    main()