import time
import os
from dist_calc import distance
from logbase2 import logbase2
# from collections import defaultdict
from collections import Counter
from math import log2
import functools
from sortedcontainers import SortedList
from operator import add

CURR_PATH = os.path.dirname(__file__)
os.chdir(CURR_PATH)

class StringUtils:
    def getOverlappingSubstrings(string: str, size: int) -> list[str]:
        substrings = [string[i:i+size] for i in range(0, len(string)-size+1)]
        return substrings

    def displayMotifs(sequences, motif, L):
        result = ""
        p = motif[0] # motif prototype
        result += f'{p.seq_idx}\t'
        result += sequences[p.seq_idx][:p.start_pos]
        result += '\033[92m' + sequences[p.seq_idx][p.start_pos:p.start_pos+L] + '\033[0m'
        result += sequences[p.seq_idx][p.start_pos+L:] + '\n'
    
        for lmer in motif[1:]:
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
        def __init__(self, lsequence, seq_idx, start_pos, w, N):
                self.lsequence: str = lsequence
                self.seq_idx = seq_idx
                self.start_pos = start_pos

                self.distances = [SortedList() for _ in range(N)]
                self.distances[self.seq_idx].add(SequenceProcessor.Instance(0, self))

                self.instance_mask = [False for _ in range(N)]
                self.instance_mask[self.seq_idx] = True
        def __str__(self) -> str:
            return f'{self.seq_idx}:[{self.lsequence}]'
        def __repr__(self) -> str:
            return self.__str__()
    
    @functools.total_ordering
    class Instance:
        def __init__(self, dist, lmer):
            self.dist = dist
            self.lmer = lmer
        def __eq__(self, other):
            return isinstance(other, type(self)) and self.dist == other.dist
        def __lt__(self, other):
            return isinstance(other, type(self)) and self.dist > other.dist
        
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
                            second_lmer.instance_mask[first_seq] = True

                            first_lmer.distances[second_seq].add(SequenceProcessor.Instance(dist, second_lmer))
                            second_lmer.distances[first_seq].add(SequenceProcessor.Instance(dist, first_lmer))
    
    def removeDegenerateLmers(lmer_list, seq_count):
        candidates = []
        for lmer in lmer_list:
            if sum(lmer.instance_mask) == seq_count:
                candidates.append(lmer)
        return candidates
    
    def IC(motif):
        lmer_strings = [l.lsequence for l in motif]
        string_length = len(lmer_strings[0])
        seq_count = len(lmer_strings)
        total = seq_count*string_length
        concat_string = ''.join(lmer_strings)
        counts = Counter(concat_string)
        A_background = (counts['A']+1)/total
        C_background = (counts['C']+1)/total
        G_background = (counts['G']+1)/total
        T_background = (counts['T']+1)/total
        information_content = 0

        for p in range(string_length):
            A_odds = sum([c=='A' for s in lmer_strings for c in s[p]])/seq_count
            C_odds = sum([c=='C' for s in lmer_strings for c in s[p]])/seq_count
            G_odds = sum([c=='G' for s in lmer_strings for c in s[p]])/seq_count
            T_odds = sum([c=='T' for s in lmer_strings for c in s[p]])/seq_count

            information_content += A_odds*logbase2(A_odds/A_background) + C_odds*logbase2(C_odds/C_background) + G_odds*logbase2(G_odds/G_background) + T_odds*logbase2(T_odds/T_background)

        return information_content
    
    def chooseBestInstance(curr_instances, possible_instances):
        # instance_strings = [s.lmer.lsequence for s in possible_instances]
        # try:
        lmer_strings = [l.lsequence for l in curr_instances]
        # except:
            # print(curr_instances)
        string_length = len(lmer_strings[0])
        seq_count = len(lmer_strings)
        total = seq_count*string_length
        concat_string = ''.join(lmer_strings)
        counts = Counter(concat_string)
        A_background = counts['A']+1
        C_background = counts['C']+1
        G_background = counts['G']+1
        T_background = counts['T']+1
        A_positional_counts = [sum([c=='A' for s in lmer_strings for c in s[p]]) for p in range(string_length)]
        C_positional_counts = [sum([c=='C' for s in lmer_strings for c in s[p]]) for p in range(string_length)]
        G_positional_counts = [sum([c=='G' for s in lmer_strings for c in s[p]]) for p in range(string_length)]
        T_positional_counts = [sum([c=='T' for s in lmer_strings for c in s[p]]) for p in range(string_length)]

        best_assignment = None
        best_score = -999
        for a in possible_instances:
            information_content = 0
            assgn_counts = Counter(a.lmer.lsequence)
            A_bg_new = (A_background + assgn_counts['A'])/(total+1)
            C_bg_new = (C_background + assgn_counts['C'])/(total+1)
            G_bg_new = (G_background + assgn_counts['G'])/(total+1)
            T_bg_new = (T_background + assgn_counts['T'])/(total+1)
            A_odds = [x/(seq_count+1) for x in list(map(add, A_positional_counts, [a.lmer.lsequence[i] == 'A' for i in range(string_length)]))]
            C_odds = [x/(seq_count+1) for x in list(map(add, C_positional_counts, [a.lmer.lsequence[i] == 'C' for i in range(string_length)]))]
            G_odds = [x/(seq_count+1) for x in list(map(add, G_positional_counts, [a.lmer.lsequence[i] == 'G' for i in range(string_length)]))]
            T_odds = [x/(seq_count+1) for x in list(map(add, T_positional_counts, [a.lmer.lsequence[i] == 'T' for i in range(string_length)]))]
            # print(A_odds)
            for p in range(string_length):
                information_content += A_odds[p]*logbase2(A_odds[p]/A_bg_new) + C_odds[p]*logbase2(C_odds[p]/C_bg_new) + G_odds[p]*logbase2(G_odds[p]/G_bg_new) + T_odds[p]*logbase2(T_odds[p]/T_bg_new)
            # if information_content < 0: print(information_content)
            if information_content > best_score:
                best_score = information_content
                best_assignment = a.lmer
        # print(best_score)
        return best_assignment
    
    def findBestMotif(candidate_lmers, seq_count):
        best_motif = None
        best_score = 0
        for cl in candidate_lmers:
            instance_list = [cl]
            seq_indices = list(range(seq_count))
            seq_indices.pop(cl.seq_idx)
            for seq_idx in seq_indices:
                if len(cl.distances[seq_idx]) == 1:
                    instance_list.append(cl)
                else:
                    instance_list.append(MotifFinder.chooseBestInstance(instance_list, cl.distances[seq_idx]))
            score = MotifFinder.IC(instance_list)
            # print(score)
            if score > best_score:
                best_score = score
                best_motif = instance_list
        return best_motif

    

def main():
    sequences = SequenceProcessor.readFromSitesFile("datasets\\MA0014.2.sites", n_seq=120, seq_len=60)
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
    motif = MotifFinder.findBestMotif(candidate_lmers, N)
    ########
    end_time = time.perf_counter_ns()

    StringUtils.displayMotifs(sequences, motif, L)
    print(f"Time to find motif: {(end_time-start_time)/(10**6)} ms")

if __name__ == "__main__":
    main()