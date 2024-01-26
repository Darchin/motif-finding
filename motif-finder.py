import time
import os
from dist_calc import distance
from logbase2 import logbase2
import numpy as np
from collections import Counter
import functools
from sortedcontainers import SortedList

CURR_PATH = os.path.dirname(__file__)
os.chdir(CURR_PATH)

class StringUtils:
    def getOverlappingSubstrings(string: str, size: int) -> list[str]:
        substrings = [string[i:i+size] for i in range(0, len(string)-size+1)]
        return substrings

    def displayMotifs(sequences, motif, L):
        instances = motif[0]
        pattern = motif[1]
        result = ""
        result += f'{pattern.seq_idx}\t'
        result += sequences[pattern.seq_idx][:pattern.start_pos]
        result += '\033[92m' + sequences[pattern.seq_idx][pattern.start_pos:pattern.start_pos+L] + '\033[0m'
        result += sequences[pattern.seq_idx][pattern.start_pos+L:] + '\n'
    
        for lmer in instances:
            result += f'{lmer.seq_idx}\t'
            result += sequences[lmer.seq_idx][:lmer.start_pos]
            result += '\033[91m' + sequences[lmer.seq_idx][lmer.start_pos:lmer.start_pos+L] + '\033[0m'
            result += sequences[lmer.seq_idx][lmer.start_pos+L:] + '\n'
        
        print(result)

class SequenceProcessor:
    def __init__(self, N, w):
        self.maximum_wiggle_room = w
        self.sequence_count = N
    class Lmer:
        def __init__(self, lsequence, seq_idx, start_pos, N):
                self.lsequence: str = lsequence
                self.seq_idx: int = seq_idx
                self.start_pos: int = start_pos

                self.distances: list[SortedList] = [SortedList() for _ in range(N)]
                self.distances[self.seq_idx].add(SequenceProcessor.Instance(0, self))

                self.instance_mask: list[bool] = [False for _ in range(N)]
                self.instance_mask[self.seq_idx] = True

                c = Counter(self.lsequence)
                self.counts: list[int] = [c['A'], c['C'], c['G'], c['T']] #A, C, G, T
                            
                self.A_mask = [char=='A' for char in lsequence]
                self.C_mask = [char=='C' for char in lsequence]
                self.G_mask = [char=='G' for char in lsequence]
                self.T_mask = [char=='T' for char in lsequence]
        def __str__(self) -> str:
            return f'{self.seq_idx}:[{self.lsequence}]'
        def __repr__(self) -> str:
            return self.__str__()
    
    @functools.total_ordering
    class Instance:
        def __init__(self, dist, lmer):
            self.dist = dist
            self.lmer = lmer
        def __contains__(self, dist):
            return dist == self.dist
        def __eq__(self, other):
            return (isinstance(other, type(self)) and self.dist == other.dist)# or (isinstance(other, int) and self.dist == other)
        def __lt__(self, other):
            return isinstance(other, type(self)) and self.dist < other.dist
        def __str__(self):
            return self.dist
        def __repr__(self):
            return f'{self.__str__()}'
        
    def createLmer(self, lsequence, seq_idx, start_pos):
        return self.Lmer(lsequence, seq_idx, start_pos, self.sequence_count)
    
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
                            second_lmer.instance_mask[first_seq] = True

                            first_lmer.distances[second_seq].add(SequenceProcessor.Instance(dist, second_lmer))
                            second_lmer.distances[first_seq].add(SequenceProcessor.Instance(dist, first_lmer))
    
    def removeDegenerateLmers(lmer_list, seq_count):
        candidates = []
        for lmer in lmer_list:
            if sum(lmer.instance_mask) == seq_count:
                candidates.append(lmer)
        return candidates
    
    def calculateIC(lmers, L):
        seq_count = len(lmers)
        bg_tot = seq_count * L

        A_bg = (sum([lmer.counts[0] for lmer in lmers])+1)/bg_tot
        C_bg = (sum([lmer.counts[1] for lmer in lmers])+1)/bg_tot
        G_bg = (sum([lmer.counts[2] for lmer in lmers])+1)/bg_tot
        T_bg = (sum([lmer.counts[3] for lmer in lmers])+1)/bg_tot

        A_pc = [sum([lmer.A_mask[p] for lmer in lmers]) for p in range(L)]
        C_pc = [sum([lmer.C_mask[p] for lmer in lmers]) for p in range(L)]
        G_pc = [sum([lmer.G_mask[p] for lmer in lmers]) for p in range(L)]
        T_pc = [sum([lmer.T_mask[p] for lmer in lmers]) for p in range(L)]

        A_odds = [(A_pc[p])/seq_count for p in range(L)]
        C_odds = [(C_pc[p])/seq_count for p in range(L)]
        G_odds = [(G_pc[p])/seq_count for p in range(L)]
        T_odds = [(T_pc[p])/seq_count for p in range(L)]
        
        information_content = 0
        for p in range(L):
            information_content += A_odds[p]*logbase2(A_odds[p]/A_bg) + C_odds[p]*logbase2(C_odds[p]/C_bg) + G_odds[p]*logbase2(G_odds[p]/G_bg) + T_odds[p]*logbase2(T_odds[p]/T_bg)

        return information_content
    
    def chooseBestInstance(curr_instances, possible_instances, L):
        # This parameter includes the new instance
        seq_count = len(curr_instances) + 1
        bg_tot = seq_count * L
        # The following parameters are excluding 'possible_instances'.
        A_curr_tot = sum([lmer.counts[0] for lmer in curr_instances])
        C_curr_tot = sum([lmer.counts[1] for lmer in curr_instances])
        G_curr_tot = sum([lmer.counts[2] for lmer in curr_instances])
        T_curr_tot = sum([lmer.counts[3] for lmer in curr_instances])

        A_curr_pc = [sum([lmer.A_mask[p] for lmer in curr_instances]) for p in range(L)]
        C_curr_pc = [sum([lmer.C_mask[p] for lmer in curr_instances]) for p in range(L)]
        G_curr_pc = [sum([lmer.G_mask[p] for lmer in curr_instances]) for p in range(L)]
        T_curr_pc = [sum([lmer.T_mask[p] for lmer in curr_instances]) for p in range(L)]

        best_assignment = None
        best_score = -999
        for instance in possible_instances:
            new_lmer = instance.lmer
            information_content = 0

            A_bg = (A_curr_tot + new_lmer.counts[0] + 1)/bg_tot
            C_bg = (C_curr_tot + new_lmer.counts[1] + 1)/bg_tot
            G_bg = (G_curr_tot + new_lmer.counts[2] + 1)/bg_tot
            T_bg = (T_curr_tot + new_lmer.counts[3] + 1)/bg_tot

            A_odds = [(A_curr_pc[p] + new_lmer.A_mask[p])/seq_count for p in range(L)]
            C_odds = [(C_curr_pc[p] + new_lmer.C_mask[p])/seq_count for p in range(L)]
            G_odds = [(G_curr_pc[p] + new_lmer.G_mask[p])/seq_count for p in range(L)]
            T_odds = [(T_curr_pc[p] + new_lmer.T_mask[p])/seq_count for p in range(L)]

            for p in range(L):
                information_content +=\
                      A_odds[p]*logbase2(A_odds[p]/A_bg) + C_odds[p]*logbase2(C_odds[p]/C_bg) +\
                      G_odds[p]*logbase2(G_odds[p]/G_bg) + T_odds[p]*logbase2(T_odds[p]/T_bg)
            if information_content > best_score:
                best_score = information_content
                best_assignment = new_lmer
        return best_assignment
    def pruneUnlikelyMotifs(lmers, N):
        sorted_avg_dist_list = SortedList([(sum([sl[0].dist for sl in lmer.distances])/N, lmer) for lmer in lmers], key=lambda x: x[0])
        # print(sorted_avg_dist_list[:20])
        mean_of_avgs = np.mean([x[0] for x in sorted_avg_dist_list])
        idx = sorted_avg_dist_list.bisect_right((mean_of_avgs, None))
        return [tup[1] for tup in sorted_avg_dist_list[:idx]]
    def pruneDistanceList(lmers, w):
        dummy_instance = SequenceProcessor.Instance(0, None)
        for lmer in lmers:
            for i in range(len(lmer.distances)):
                closest = lmer.distances[i][0].dist
                dummy_instance.dist = closest + w + 1
                if dummy_instance in lmer.distances[i]:
                    drop_idx = lmer.distances[i].index(dummy_instance)
                else: continue
                del lmer.distances[i][drop_idx:]
    def findBestMotif(candidate_lmers, seq_count, L):
        best_motif = None
        best_score = 0
        for cl in candidate_lmers:
            instance_list = []
            for seq_idx in range(seq_count):
                if seq_idx == cl.seq_idx: instance_list.append(cl)
                elif len(cl.distances[seq_idx]) == 1: instance_list.append(cl.distances[seq_idx][0].lmer)
                else: 
                    best_inst = MotifFinder.chooseBestInstance(instance_list, cl.distances[seq_idx], L)
                    instance_list.append(best_inst)

            score = MotifFinder.calculateIC(instance_list, L)
            if score > best_score:
                best_score = score
                best_motif = (instance_list, cl)
        return best_motif

def main():
    sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0014.2.sites", n_seq=60, seq_len=60)
    L = 10
    d = 5
    w = 1
    N = len(sequences)
    sq = SequenceProcessor(N, w)
    lmer_seq_lst = sq.extractLmersFromSequenceList(sequences, L)
    lmer_flatlist = [lmer for lmer_list in lmer_seq_lst for lmer in lmer_list]
    
    start_time = time.perf_counter_ns()
    ########
    MotifFinder.findInstances(lmer_seq_lst, N, d)
    candidate_lmers = MotifFinder.removeDegenerateLmers(lmer_flatlist, N)
    MotifFinder.pruneDistanceList(candidate_lmers, w)
    likely_patterns = MotifFinder.pruneUnlikelyMotifs(candidate_lmers, N)
    motif = MotifFinder.findBestMotif(likely_patterns, N, L)
    ########
    end_time = time.perf_counter_ns()

    StringUtils.displayMotifs(sequences, motif, L)
    correct_instances = sum([(inst.start_pos == 50) for inst in motif[0]])
    accuracy = correct_instances/N

    # for MA0014.2.sites
    print(f'Time to find motif: {(end_time-start_time)/(10**6)} ms')
    print(f'Instances correctly Identified: {correct_instances}/{N}')
    print(f'Accuracy: {accuracy*100}%')

if __name__ == "__main__":
    main()